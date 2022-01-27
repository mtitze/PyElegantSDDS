# -*- coding: utf-8 -*-

"""
Module pyelegantsdds.elegantrun 
=================================================================

A module containing the class ElegantRun to run Elegant simulations in 
a singularity container.

"""
from .elegant_command import ElegantCommandFile
from .tools.shell_process import call

# ==============================================================================
#
# HELPER FUNCTIONS
#
# ==============================================================================


def write_parallel_elegant_script(rootname='temp'):
    """
    Method to generate a script that runs
    pelegant from bash.
    """

    bashstr = f'''#!/usr/bin/env bash
if [ $# == 0 ] ; then
    echo "usage: {rootname}_sif_pelegant <inputfile>"
    exit 1
fi
n_cores=`grep processor /proc/cpuinfo | wc -l`
echo The system has $n_cores cores.
n_proc=$((n_cores-1))
echo $n_proc processes will be started.
if [ ! -e ~/.mpd.conf ]; then
    echo "MPD_SECRETWORD=secretword" > ~/.mpd.conf
    chmod 600 ~/.mpd.conf
fi
mpiexec -host $HOSTNAME -n $n_proc Pelegant  $1 $2 $3 $4 $5 $6 $7 $8 $9
'''

    # write to file
    with open(f"{rootname}_sif_pelegant.sh", "w") as f:
        f.write(bashstr)


def write_parallel_run_script(sif, rootname='temp'):
    """
    Method to generate parallel elegant run
    script.

    Parameters:
    -----------
    sif: str path to singularity container
    """
    
    bashstr = f'''#!/bin/bash
pele={sif}
cmd="bash {rootname}_sif_pelegant.sh"
        
$pele $cmd $1
'''

    # write to file
    with open(f"{rootname}_run_pelegant.sh", "w") as f:
        f.write(bashstr)

# ==============================================================================
#
# MAIN CLASS
#
# ==============================================================================


class ElegantRun:
    """
    Class to interact with Elegant and Parallel Elegant from Python.
    """

    def __init__(self, sif, lattice: str, energy_gev: float, beamline: str, parallel=False, rootname='run'):
        self.sif = sif
        self.lattice = lattice
        self.energy_gev = energy_gev
        self.beamline = beamline
        self.parallel = parallel
        
        self.set_rootname(rootname)
        self.commandfile = ElegantCommandFile(self.commandfile_name)
        
        # setting up executable
        if parallel:
            self._write_parallel_script()
            self.exec = "bash {}".format(self.pelegant)
        else:
            self.exec = "{} elegant ".format(self.sif)            
            
    def set_rootname(self, rootname):
        """
        This function initializes the commandfile and can be used to change the rootname conveniently.
        """
        self.rootname = rootname
        self.commandfile_name = f"{rootname}.ele"

    def _write_parallel_script(self):
        """
        Method sets self.pelegant to the script file.
        """
        write_parallel_elegant_script(rootname=self.rootname)
        write_parallel_run_script(self.sif, rootname=self.rootname)
        self.pelegant = f"{self.rootname}_run_pelegant.sh"
        
    def copy(self, **kwargs):
        '''
        Return a copy of this instance (without command history)
        '''
        return self.__class__(sif=kwargs.get('sif', self.sif), 
                              lattice=kwargs.get('lattice', self.lattice),
                              energy_gev=kwargs.get('energy_gev', self.energy_gev),
                              beamline=kwargs.get('beamline', self.beamline),
                              parallel=kwargs.get('parallel', self.parallel),
                              rootname=kwargs.get('rootname', self.rootname))

    def clearCommands(self):
        """Clear the command list."""
        self.commandfile.clear()

    def clearCommandHistory(self):
        """Clear the command history."""
        self.commandfile.clearHistory()

    def clearAll(self):
        """Clears both command list and command history."""
        self.clearCommands()
        self.clearCommandHistory()

    def run(self, verbose=True, **kwargs):
        """
        Run the commandfile.
        """
        # check if commandfile is not empty
        if len(self.commandfile.commandlist) == 0:
            print("Commandfile empty - nothing to do.")
            return

        # write Elegant command file to disk
        self.commandfile.write()
        self.commandfile.clear()

        cmdstr = "{} {}.ele".format(self.exec, self.rootname)
        call(cmdstr, rootname=self.rootname, verbose=verbose)
            
    ######################     
    # Elegant commands
    ######################
    # The following commands shorten the existing elegant commands.
            
    def add_run_setup(self, **kwargs):
        """
        Add run_setup command.
        """
        self.commandfile.addCommand(
            "run_setup",
            lattice=self.lattice,
            use_beamline=self.beamline,
            p_central_mev=self.energy_gev*1e3,
            # centroid="%s.cen",
            default_order=kwargs.get("default_order", 2),
            concat_order=kwargs.get("concat_order", 0),
            rootname=self.rootname,
            losses="%s.lost",
            losses_include_global_coordinates=kwargs.get("losses_include_global_coordinates", 0),           
            #losses_s_limit=f'{kwargs.get("losses_s_limit_1", -9999)} {kwargs.get("losses_s_limit_2", 9999)}',
            acceptance="%s.acc",
            parameters="%s.params",
            semaphore_file="%s.done",
            magnets="%s.mag",  # for plotting profile
            final="%s.fin",
            output="%s.out",
        )

    def add_twiss_output(self, **kwargs):
        """
        Add basic twiss.
        
        Resonance driving terms (according to https://ops.aps.anl.gov/manuals/elegant_latest/elegantsu79.html#x87-860007.69):
        Term Name	Explanation
        h11001 	drives x chromaticity
        h00111 	drives y chromaticity
        h20001 	drives synchro-betatron resonances
        h00201 	drives momentum-dependence of beta functions
        h10002 	drives second order dispersion
        h21000 	drives νx
        h30000 	drives 3νx
        h10110 	drives νx
        h10020 	drives νx - 2νy
        h10200 	drives νx + 2νy
        h22000 	drives dνx∕dJx
        h11110 	drives dνx∕dJy
        h00220 	drives dνy∕dJy
        h31000 	drives 2νx
        h40000 	drives 4νx
        h20110 	drives 2νx
        h11200 	drives 2νy
        h20020 	drives 2νx - 2νy
        h20200 	drives 2νx + 2νy
        h00310 	drives 2νy
        h00400 	drives 4νy
        """
        concat_order = kwargs.get('concat_order', 3)
        higher_order_chromaticity = kwargs.get('higher_order_chromaticity', 1)
        if concat_order < 2 and higher_order_chromaticity != 0:
            raise RuntimeError('If higher_order_chromaticity != 0, then it is necessary to set concat_order >= 3.') # see https://ops.aps.anl.gov/manuals/elegant_latest/elegantsu79.html#x87-860007.69
        
        self.commandfile.exists('run_setup')
        self.commandfile.addCommand(
            "twiss_output", filename=kwargs.get('filename', '%s.twi'), 
            matched=kwargs.get('matched', 1),
            radiation_integrals=kwargs.get('radiation_integrals', 1),
            concat_order=concat_order,
            higher_order_chromaticity=higher_order_chromaticity,
            compute_driving_terms=kwargs.get('compute_driving_terms', 1)
        )

    def add_closed_orbit(self, **kwargs):
        '''
        Add closed_orbit command.
        '''
        self.commandfile.addCommand(
            "closed_orbit",
            output=kwargs.get('output', '%s.clo'),
            tracking_turns=kwargs.get('tracking_turns', 256)
        )
        
    def add_moments_output(self, **kwargs):
        '''
        Add moments_output command.
        '''
        # N.B:
        # parameter tracking_based_diffusion_matrix_particles may require update of Elegant ...

        self.commandfile.addCommand(
            "moments_output",
            filename=kwargs.get('filename', '%s.mom'),
            output_at_each_step=kwargs.get('output_at_each_step', 0),
            output_before_tune_correction=kwargs.get('output_before_tune_correction', 0),
            final_values_only=kwargs.get('final_values_only', 0),
            matched=kwargs.get('matched', 1),
            equilibrium=kwargs.get('equilibrium', 1),
            radiation=kwargs.get('radiation', 1),
            n_slices=kwargs.get('n_slices', 10)
        )
        
    def add_vary_element(self, **kwargs):
        """
        Add single vary element line.
        """
        self.commandfile.addCommand(
            "vary_element",
            name=kwargs.get("name", "*"),
            item=kwargs.get("item", "L"),
            intial=kwargs.get("initial", 0.0000),
            final=kwargs.get("final", 0.0000),
            index_number=kwargs.get("index_number", 0),
            index_limit=kwargs.get("index_limit", 1),
            start_occurence=kwargs.get('start_occurence', -1),
            end_occurence=kwargs.get('end_occurence', -1)
        )
        
    def add_alter_elements(self, **kwargs):
        """Add alter_element command."""
        self.commandfile.exists('run_setup')
        self.commandfile.addCommand("alter_elements", **kwargs)
        
    def add_rf_setup(self, **kwargs):
        """Add rf_setup command."""
        self.commandfile.exists('twiss_output')
        self.commandfile.addCommand(
            "rf_setup",
            filename="%s.rf", 
            harmonic=kwargs.get('harmonic', 0),
            total_voltage=kwargs.get('total_voltage', 0),
            near_frequency=kwargs.get('near_frequency', 0),
            track_for_frequency=kwargs.get('track_for_frequency', 0)     
        )
        
    def add_momentum_aperture(self, **kwargs):
        self.commandfile.addCommand(
            "momentum_aperture",
            output=kwargs.get('output', '%s.mmap'),
            x_initial=kwargs.get('x_initial', 0),
            y_initial=kwargs.get('y_initial', 0),  # The initial x and y coordinate values for tracking. It is essential that y_initial be nonzero if one wants to see losses due to vertical resonances. 
            fiducialize=kwargs.get('fiducialize', 0), # If given, an initially on-energy particle is tracked before the momentum aperture search begins, in order to fiducialize the reference momentum. This is useful if there are synchrotron radiation losses or energy gain due to cavities in the system. 
            delta_negative_start = kwargs.get('delta_negative_start', 0.0),
            delta_positive_start = kwargs.get('delta_positive_start', 0.0),
            delta_negative_limit = kwargs.get('delta_negative_limit', -0.1),
            delta_positive_limit = kwargs.get('delta_positive_limit', 0.1),
            delta_step_size = kwargs.get('delta_step_size', 0.01),
            include_name_pattern = kwargs.get('include_name_pattern', None),
            include_type_pattern = kwargs.get('include_type_pattern', None),
            output_mode = kwargs.get('output_mode', 0) # Normally, elegant puts the values for positive and negative momentum aperture in different columns. Each element thus has a single row of data in the output file. If output_mode=1, elegant instead puts the values for positive and negative apertures in successive rows, with a reduced number of columns. This is mostly advantageous for the parallel version, since it allows using twice as many simultaneous processors. If output_mode=2, elegant tracks many more probe particles simultaneously, which is better for massively parallel systems. The number of particles tracked is the number of elements selected times the number of probe points between delta_negative_limit and delta_positive_limit. 
        )
        
    def add_tscatter_elements(self, **kwargs):
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
        
        # N.B. add_at_end = 1 and add_at_start = 1 are required due to this update:
        # https://ops.aps.anl.gov/manuals/elegant_latest/elegantsu260.html#x274-42800013.23
        
        
        self.commandfile.addCommand(
            "insert_elements",
            name=kwargs.get("name", "*"),
            type=kwargs.get("type", "*[DBQS]*"),
            exclude="",
            s_start=kwargs.get("s_start", -1),
            s_end=kwargs.get("s_end", -1),
            skip=kwargs.get("skip", 1),
            insert_before=kwargs.get("insert_before", 0),
            add_at_end=kwargs.get("add_at_end", 1),
            add_at_start=kwargs.get("add_at_start", 1),
            element_def=kwargs.get("element_def", r'"TSC: TSCATTER"'))

    def add_touschek_scatter(self, **kwargs):
        '''
        Add touschek_scatter command to the command list.
        '''
        
        self.commandfile.addCommand(
            "touschek_scatter",
            charge=kwargs.get('charge', 0),
            frequency=kwargs.get('frequency', 1),
            emit_x=kwargs.get('emit_x', 0),
            #emit_nx=kwargs.get('emit_nx', 0),
            emit_y=kwargs.get('emit_y', 0),
            #emit_ny=kwargs.get('emit_ny', 0),
            sigma_dp=kwargs.get('sigma_dp', 0),
            sigma_s=kwargs.get('sigma_s', 0),
            Momentum_Aperture_scale=kwargs.get('Momentum_Aperture_scale', 0.85),
            Momentum_Aperture=kwargs.get('Momentum_Aperture', None),
            XDist=kwargs.get('XDist', None),
            YDist=kwargs.get('YDist', None),
            ZDist=kwargs.get('ZDist', None),
            TranDist=kwargs.get('TranDist', None),
            FullDist=kwargs.get('FullDist', None),
            bunch=kwargs.get('bunch', '%s-%03ld.bun'),
            loss=kwargs.get('loss', '%s-%03ld.los'),
            distribution=kwargs.get('distribution', '%s-%03ld.dis'),
            initial=kwargs.get('initial', '%s-%03ld.ini'),
            output=kwargs.get('output', '%s-%03ld.out'),
            nbins=kwargs.get('nbins', 100),
            sbin_step=kwargs.get('sbin_step', 1),
            n_simulated=kwargs.get('n_simulated', 5000000),
            ignored_portion=kwargs.get('ignored_portion', 0.01),
            i_start=kwargs.get('i_start', 0),
            i_end=kwargs.get('i_end', 1),
            do_track=kwargs.get('do_track', 0),
            match_position_only=kwargs.get('match_position_only', 0),
            overwrite_files=kwargs.get('overwrite_files', 1),
            verbosity=kwargs.get('verbosity', 0)
        )
        
    def add_floor_coordinates(self, **kwargs):
        '''
        Add floor_coordinates command.
        '''
        self.commandfile.addCommand(
            "floor_coordinates",
            filename='%s.flr',
            X0=kwargs.get('x0', 0.0),  
            Z0=kwargs.get('z0', 0.0),  
            theta0=kwargs.get('theta0', 0.0),
            include_vertices=kwargs.get('include_vertices', 0),  
            vertices_only=kwargs.get('vertices_only', 0),  
            magnet_centers=kwargs.get('magnet_centers', 0), 
            store_vertices=kwargs.get('store_vertices', 0)  
        )
        
    def add_frequency_map(self, **kwargs):
        """
        Add elegant standard fma command.
        """

        self.commandfile.addCommand(
            "frequency_map",
            output="%s.fma",
            xmin=kwargs.get("xmin", -0.1),
            xmax=kwargs.get("xmax", 0.1),
            ymin=kwargs.get("ymin", 1e-6),
            ymax=kwargs.get("ymax", 0.1),
            delta_min=kwargs.get("delta_min", 0),
            delta_max=kwargs.get("delta_max", 0),
            nx=kwargs.get("nx", 21),
            ny=kwargs.get("ny", 21),
            ndelta=kwargs.get("ndelta", 1),
            verbosity=0,
            include_changes=kwargs.get("include_changes", 1),
            quadratic_spacing=kwargs.get("quadratic_spacing", 0),
            full_grid_output=kwargs.get("full_grid_output", 1)
        )