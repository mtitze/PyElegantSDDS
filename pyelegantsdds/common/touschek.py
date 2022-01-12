from ..tools.shell_process import call
from ..sdds import SDDS

from .main import ElegantRunToolkit
from .observe import watchpoints
from .rf import rftools
from .radiation import radiation

from shutil import copyfile

class touschek(ElegantRunToolkit):
    
    def piwinski(self, outstr='tlife', verbose=True, **kwargs):
        '''
        Compute the Touschek lifetime parameters from an Elegant dynamic momentum map.
        The parameters are found in Piwinski's analytical formula.
        
        Note: If using Pelegant to compute the momentum aperture with output_mode=1, 
              it is necessary to first run the script reorganizeMmap to put the data into the form needed by touschekLifetime.
              
        Parameters
        ----------
        outstr: str
            The file name appendix of the SDDS file containing the Touschek lifetime data.
            
        **kwargs
            Arguments passed to touschekLifetime routine 
            (see https://ops.aps.anl.gov/manuals/elegant_latest/elegantsu105.html#x114-1130008.23)
        '''
        # Own note regarding the output SDDS:
        # I guess 'P' and 'N' refer to 'positive' and 'negative' values according to what is done in the
        # calculation of the dynamic momentum aperture. Accordingly, TmP and TnP are related to the local
        # life times, and FP, FN belong to the function F in Piwinski's paper. The entries B1 and B2 belong to
        # the respective functions in Piwinski's paper.
        
        # set cmdstr and run
        cmdstr = f'{self.er.sif} touschekLifetime {self.er.rootname}.{outstr}'
        
        # put in some default values in case the user did not provided them (for convenience)
        if 'twiss' not in kwargs.keys():
            kwargs['twiss'] = f'{self.er.rootname}.twi'
        if 'aper' not in kwargs.keys():
            kwargs['aper'] = f'{self.er.rootname}.mmap'

        # add the remaining keys to the command string
        for key in kwargs.keys():
            cmdstr += f' -{key}={kwargs[key]}'
        
        # execute the command
        call(cmdstr, rootname=self.er.rootname, verbose=verbose)
        
        # load output from the SDDS Touschek-lifetime file
        return SDDS(self.er.sif, f"{self.er.rootname}.{outstr}", 0)
        
    
    def add_tscatter(self, **kwargs):
        """Add TSCATTER element(s) into the lattice.
        
        By default, a TSCATTER element is inserted in the lattice after each drift, bend, quadrupole and
        sextupole element.
        
        Further details here:
        https://www3.aps.anl.gov/forums/elegant/viewtopic.php?f=14&t=476&sid=467516201e9bd336985a828b20127b0b
        
        Parameters
        ----------
        save: str, optional
            If given, lattice will be stored after the insertion of the TSCATTER elements.
        """
        self.er.commandfile.addCommand(
            "insert_elements",
            name=kwargs.get("name", "*"),
            type=kwargs.get("type", "*[DBQS]*"),
            exclude="",
            s_start=kwargs.get("s_start", -1),
            s_end=kwargs.get("s_end", -1),
            skip=kwargs.get("skip", 1),
            insert_before=kwargs.get("insert_before", 0),
            add_at_end=kwargs.get("add_at_end", 0),
            add_at_start=kwargs.get("add_at_start", 1),
            element_def=kwargs.get("element_def", r'"TSC: TSCATTER"'))
        
        
    def prepare_tscatter(self, n_passes, add_watch_start=True, rf=0, rad=True, verbose=True, **kwargs):
        '''
        Prepare the given lattice to be used with TSCATTER elements: TSCATTER elements are inserted into
        the lattice and the momentum aperture will be calculated.
        
        Parameters
        ----------
        n_passes: int
            Number of revolutions to be tracked while determining the dynamic momentum aperture. This
            should be sufficiently large to cover at least one synchrotron period.
        '''
        
        # construct command file
        self.er.commandfile.clear()
        self.er.add_run_setup(**kwargs)
        
        # tmp, may extract from lattice using rf_handling
        assert 'total_voltage' in kwargs.keys()
        assert 'rf_frequency' in kwargs.keys()
        
        if rf != 0:
            rfc = rftools(self.er)
            rfc.rf_handling(rf_status=rf, **kwargs)
            
        if rad:
            radc = radiation(self.er)
            radc.add_radiation_damping(**kwargs)
        
        if add_watch_start:
            msc = watchpoints(self.er)
            msc.add_watch_at_start()

        # add elements at which to compute the Touschek scattering events
        self.add_tscatter(**kwargs)
        
        # save modified lattice to new file for later use in self.run_tscatter
        self.tsc_lattice_name = f'{self.er.rootname}_tsc.lte'    
        self.er.commandfile.addCommand("save_lattice",
                                        filename=self.tsc_lattice_name,
                                        output_seq=1)
        
        # perform a twiss calculation
        self.tsc_twiss_name = f'{self.er.rootname}_tsc.twi'
        self.er.add_twiss_output(matched=1, radiation_integrals=1,
                                 filename=self.tsc_twiss_name)
            
        # compute the dynamic momentum aperture
        self.er.commandfile.addCommand("run_control", n_passes=n_passes)
        self.er.add_momentum_aperture(**kwargs)
            
        # run the above commands in Elegant
        self.er.run()
        
        # copy the resulting data containing the dynamic momentum aperture to a different to a specific name in order to prevent that it may get overwritten by some other routine (e.g. the dynap routine).
        self.tsc_dynmomap_file = f"{self.er.rootname}_tsc.mmap"
        copyfile(f"{self.er.rootname}.mmap", self.tsc_dynmomap_file)
        
        # load the SDDS aperture data & sort them according to s
        sdds_mmap = SDDS(self.er.sif, self.tsc_dynmomap_file, 0, rootname=self.er.rootname)
        sdds_mmap.sort() # sort by s (default). Required if running Pelegant (see https://ops.aps.anl.gov/manuals/elegant_latest/elegantsu50.html#x58-570007.40)
        # store the results locally for convenience
 
        self.tsc_dynmomap_data = sdds_mmap.getColumnValues()
        
        # load the SDDS twiss data & sort them according to s
        sdds_twiss = SDDS(self.er.sif, self.tsc_twiss_name, 0, rootname=self.er.rootname)
        sdds_twiss.sort() # sort by s (default).
        # store the results locally for convenience
        self.tsc_twiss = sdds_twiss.getColumnValues()
        
        
    def run_tscatter(self, n_passes, verbose=True, **kwargs):
        
        if not 'tsc_lattice_name' in kwargs.keys():
            assert hasattr(self, 'tsc_lattice_name'), "A lattice file needs to be provided. Run self.prepare_tscatter or provide the file by the 'tsc_lattice_name' argument."
            tsc_lattice_name = self.tsc_lattice_name
        else:
            tsc_lattice_name = kwargs['tsc_lattice_name']
            
        if verbose:
            print (f'Touschek scatter run using lattice:\n {tsc_lattice_name}')
            
        
        
        '''self.er.commandfile.addCommand("run_control", n_passes=n_passes)

        #self.commandfile.addCommand("bunched_beam")
        self.er.commandfile.addCommand(
            "sdds_beam",
            input=sdds_beam_file,
            input_type='"elegant"',
        )'''
            
