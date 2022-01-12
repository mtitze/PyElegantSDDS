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
        self.er.add_tscatter_elements(**kwargs)
        
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
        
        
    def run_tscatter(self, n_passes, sdds_beam_file, verbose=True, **kwargs):
        
        # check input consistency
        assert 'charge' in kwargs.keys(), 'charge must be provided.'
        assert kwargs['charge'] != 0, 'charge must be non-zero.'
        
        assert 'emit_x' in kwargs.keys() or 'emit_nx' in kwargs.keys(), 'emit_x or emit_nx must be provided.'
        assert 'emit_y' in kwargs.keys() or 'emit_ny' in kwargs.keys(), 'emit_y or emit_ny must be provided.'
        
        assert 'sigma_s' in kwargs.keys(), 'sigma_s must be provided.'
        assert 'sigma_dp' in kwargs.keys(), 'sigma_dp must be provided.'
        
        if not 'tsc_lattice_name' in kwargs.keys():
            assert hasattr(self, 'tsc_lattice_name'), "A lattice file needs to be provided. Run self.prepare_tscatter or provide the file by a 'tsc_lattice_name' argument."
            tsc_lattice_name = self.tsc_lattice_name
        else:
            tsc_lattice_name = kwargs['tsc_lattice_name']
            
        if not 'Momentum_Aperture' in kwargs.keys():
            assert hasattr(self, 'tsc_dynmomap_file'), "A dynamic momentum aperture file needs to be provided. Run self.prepare_tscatter or provide the file by a 'Momentum_Aperture' argument."
            kwargs['Momentum_Aperture'] = self.tsc_dynmomap_file
            
        if verbose:
            print (f'Touschek scatter run using lattice:\n {tsc_lattice_name}')
            
        # copy current ElegantRun parameters, now using the lattice containing the TSC elements 
        er_for_tsc = self.er.copy(lattice=tsc_lattice_name)
       
        er_for_tsc.add_run_setup(**kwargs) # N.B. only 3 keys are checked here at the moment, and they do not agree with any keys in the add_tscatter command below.
       
        # a twiss calculation is necessary prior to performing touschek_scatter
        er_for_tsc.add_twiss_output(matched=1, radiation_integrals=1)
      
        # define input for scatter simulation
        er_for_tsc.commandfile.addCommand("run_control", n_passes=n_passes)
        #self.commandfile.addCommand("bunched_beam")
        er_for_tsc.commandfile.addCommand(
            "sdds_beam",
            input=sdds_beam_file,
            input_type='"elegant"',
        )

        # add the touschek_scatter command
        er_for_tsc.add_touschek_scatter(**kwargs)

        # run everything
        er_for_tsc.run()


