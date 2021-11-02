from .elegantrun import ElegantRun
from .sdds import SDDS, SDDSCommand
import numpy as np
from scipy import constants as const
import pandas as pd
import os
import copy


class ElegantRunToolkit:
    '''
    Class to perform advanced operations using the ElegantRun functionality.
    '''
    def __init__(self, er: ElegantRun):
        self.er = er

    def create_new_er(self, new_rootname=''):
        # create new instance of current elegant run, using the same init parameters
        # as the original one, with a new rootname.
        if len(new_rootname) == 0:
            new_rootname = self.er.rootname + '_copy'
        return ElegantRun(sif=self.er.sif, lattice=self.er.lattice, energy_gev=self.er.energy_gev,
                             beamline=self.er.beamline, rootname=new_rootname,
                             parallel=self.er.parallel)
        
    def create_copy(self, new_rootname=''):
        # create a copy of the current elegant run, including its command list
        er2 = self.create_new_er(new_rootname)
        er2.commandfile.commandlist = copy.deepcopy(self.er.commandfile.commandlist)
        er2.commandfile.history = copy.deepcopy(self.er.commandfile.history)
        return er2

class rftools(ElegantRunToolkit):
    '''
    Class to perform RF modifications to an Elegant run.
    '''
        
    def add_rf_cavity(self, length, voltage, frequency, phase, **kwargs):
        """Add generic rf cavity (at lattice start by default)."""
        self.er.commandfile.exists('run_setup')
        self.er.commandfile.addCommand(
            "insert_elements",
            name=kwargs.get("name", "RF"),
            type=kwargs.get("type", "RFCA"),
            exclude="",
            s_start=kwargs.get("s_start", -1),
            s_end=kwargs.get("s_end", -1),
            skip=kwargs.get("skip", 1),
            insert_before=kwargs.get("insert_before", 0),
            add_at_end=kwargs.get("add_at_end", 0),
            add_at_start=kwargs.get("add_at_start", 1),
            #element_def=kwargs.get("element_def", r'"RF: RFCA, "')
            element_def=f'"RF: RFCA, L={length}, VOLT={voltage}, FREQ={frequency}, PHASE={phase}"'
            )

    def get_rf_get_freq_and_phase(self, total_voltage: float, harmonic: int, rad: bool=False, **kwargs):
        """Add simulation single turn to get synchronuous freq and phase to be used in next simulation.
        Produces temp.twi and temp.out that can be used with tags.
        Parameters
        ----------
        total_voltage : float, optional
            Total sum of RF voltages, by default 4*375e3
        harmonic : int, optional
            harmonic number, by default 400
        rad : bool, optional
            radiation flag, by default False
        """
        
        # create a copy of the current Elegant run to get the RF frequency and the phase.
        rootname2 = f'{self.er.rootname}_rf_freq_and_phase'
        er2 = self.create_new_er(rootname2)
        er2.add_run_setup(**kwargs)
        if rad:
            rad_er2_tools = radiation(er2)
            rad_er2_tools.add_radiation_damping(**kwargs)

        er2.add_twiss_output()
        er2.add_rf_setup(total_voltage=total_voltage, harmonic=harmonic, track_for_frequency=1)
        er2.commandfile.addCommand("run_control")
        er2.commandfile.addCommand("track")
        er2.run()
        return rootname2
        
    def rf_handling(self, rf_status, total_voltage, one_turn_rootname='', **kwargs):
        
        # rf_status: 0 no RF changes
        # rf_status: 1: insert generic RF cavity
        # rf_status: 2: use/modify existing RF cavitie(s)
        
        if rf_status == 1:
            harmonic = kwargs.get('harmonic')
            rf_length = kwargs.get('rf_length', 0)
            rf_freq = kwargs.get('rf_frequency')
            assert rf_freq != None
            rf_phase = kwargs.get('rf_phase', 180)
            
            if kwargs.get("rad", False) and len(one_turn_rootname) > 0:
                twissfile = f"{one_turn_rootname}.twi"
                print (f'Loading U0 from {twissfile} to modify rf phase ...')
                sdds_one_turn_twi = SDDS(self.er.sif, twissfile, 0, rootname=f'{one_turn_rootname}')
                twiss_params = sdds_one_turn_twi.getParameterValues()
                U0 = twiss_params.U0*1e6 # convert MeV to eV
                phase_correction = -180*np.arcsin(U0/total_voltage)
            else:
                phase_correction = 0
            rf_phase += phase_correction
            print (f'RF phase: {rf_phase} deg')
            
            self.er.add_twiss_output()
            self.add_rf_cavity(length=rf_length, voltage=total_voltage, frequency=rf_freq, phase=rf_phase, **kwargs)
            #self.add_rf_setup(filename="%s.rf", **kwargs)
            #self.add_alter_elements(name="*", type="RFCA", item="PHASE", value=f'"({phase})"') # bugfix: need to change PHASE sepeartely, because rf_setup apparently computes the frequency not independently!
        
        if rf_status == 2:
            harmonic = kwargs.get('harmonic')
            n_cavities = kwargs.get('n_cavities')
            assert n_cavities > 0 # may need to get this from the lattice automatically...
            assert harmonic != None

            # get frequency
            if len(one_turn_rootname) == 0:
                print ('one_turn_rootname not provided. Attempting to obtain the frequency from tracking ...')
                one_turn_rootname = self.get_rf_get_freq_and_phase(total_voltage=total_voltage, harmonic=harmonic, rad=False) # run with track_for_frequency=1; TODO: need to discuss case rad=True
            
            self.er.add_twiss_output()
            self.er.commandfile.addCommand("rpn_load", filename=f'{one_turn_rootname}.out', tag="ref", use_row=0) ### filename exists???

            # get energy loss
            self.er.commandfile.addCommand("rpn_load", filename=f"{one_turn_rootname}.twi", tag="reftwi", load_parameters=1)

            # update RF settings
            # VOLTAGE
            self.er.add_alter_elements(name="*", type="RFCA", item="VOLT", value=total_voltage/n_cavities)
            # FREQUENCY
            # N.B. Elegant uses RPN (reverse polnic notation) for formula expressions.
            # The definition of the operator 'rec' (= reciprocal) is hereby given in the file
            # mti/containers/elegant/3rdparty/defns.rpn
            # The following operation sets the frequencies in all RF-cavities according to 1/t_ref*harmonic.
            # !!! Note that in reality, these cavities may operate with different frequencies.
            self.er.add_alter_elements(name="*", type="RFCA", item="FREQ", value=f'"(ref.t rec {harmonic} * )"') # bugfix: brackets and " were in the wrong order!

            # PHASE
            # Note: reftwi.U0 is the energy loss per turn, in Elegant given in MeV
            self.er.add_alter_elements(name="*", type="RFCA", item="PHASE", value=f'"(180 reftwi.U0 1e6 * {total_voltage} / dasin -)"') # N.B. dasin: degree-arcus-sinus: dasin() = 180/pi*arcsin() yields degreees
            # PHASE 2nd varation
            # self.add_alter_elements(name="*", type="RFCA", item="PHASE", value='"(180)"')

            #self.add_twiss_output()
    

        
class radiation(ElegantRunToolkit):

    def add_radiation_damping(self, isr=1, synch_rad=1, isr1part=1, use_rad_dist=1, **kwargs):
        """Add radiation damping for CSBEND, KQUAD and KSEXT.
        
        isr: also set ISR (include incoherent synchrotron radiation) to the given value.
        """
        self.er.add_alter_elements(name="*", type="CSBEND", item="SYNCH_RAD", value=synch_rad)
        self.er.add_alter_elements(name="*", type="CSBEND", item="USE_RAD_DIST", value=use_rad_dist)
        self.er.add_alter_elements(name="*", type="CSBEND", item="ISR", value=isr)
        self.er.add_alter_elements(name="*", type="KQUAD", item="SYNCH_RAD", value=synch_rad)
        self.er.add_alter_elements(name="*", type="KQUAD", item="ISR", value=isr)
        self.er.add_alter_elements(name="*", type="KQUAD", item="ISR1PART", value=isr1part)
        self.er.add_alter_elements(name="*", type="KSEXT", item="SYNCH_RAD", value=synch_rad)
        self.er.add_alter_elements(name="*", type="KSEXT", item="ISR", value=isr)
        self.er.add_alter_elements(name="*", type="KSEXT", item="ISR1PART", value=isr1part)


class tracking(ElegantRunToolkit):
    
    def __init__(self, er, sdds_beam_file, **kwargs):
        ElegantRunToolkit.__init__(self, er=er, **kwargs)
        self.sdds_beam_file = sdds_beam_file # tracking requires an SDDS beam input file
    
    def add_basic_controls(self):
        """Adding basic controls for tracking"""
        # add controls
        self.er.commandfile.addCommand("run_control")
        self.er.commandfile.addCommand("bunched_beam")
        self.er.commandfile.addCommand("track")
    
    def track(self, n_passes, add_watch_start=True, rad=False, rf=0, **kwargs):
        """
        Track a set of particles.
        """
        # construct command file
        self.er.commandfile.clear()
        self.er.add_run_setup(**kwargs)
        if add_watch_start:
            msc = watchpoints(self.er)
            msc.add_watch_at_start()
        self.er.add_twiss_output()
        
        if rf != 0:
            rfc = rftools(self.er)
            rfc.rf_handling(rf_status=rf, **kwargs)
            
        if rad:
            radc = radiation(self.er)
            radc.add_radiation_damping(**kwargs)

        self.er.commandfile.addCommand("run_control", n_passes=n_passes)
            
        if kwargs.get("momap", False):
            self.er.add_momentum_aperture(**kwargs)

        #self.commandfile.addCommand("bunched_beam")
        self.er.commandfile.addCommand(
            "sdds_beam",
            input=self.sdds_beam_file,
            input_type='"elegant"',
        )

        self.er.commandfile.addCommand("track")
        
        
    #################### old code, not yet updated to new version

    def simple_single_particle_track(self, coord=np.zeros((5, 1)), **kwargs):
        """
        Track a single particle with given initial coordinates.

        Important:
        ----------
        Be careful with giving the 6th coordinate, this is beta * gamma. If not
        given it will be calculated automatically either using standard 1700 MeV
        or kwargs["pcentralmev"].
        """
        # generate particle input file
        pg = particle_generation(er=self.er)
        sdds_beam_file = pg.generate_sdds_particle_inputfile(
            man_ranges={k: v for k, v in zip(range(coord.shape[0] + 1), coord)}, **kwargs
        )

        # construct command file
        self.er.commandfile.clear()
        self.er.add_run_setup(**kwargs)
        self.er.add_watch_at_start()

        self.er.commandfile.addCommand("run_control", n_passes=kwargs.get("n_passes", 2 ** 8))
        #self.commandfile.addCommand("bunched_beam")
        self.er.commandfile.addCommand(
            "sdds_beam",
            input=sdds_beam_file,
            input_type='"elegant"',
        )
        self.er.commandfile.addCommand("track")

        # run will write command file and execute it
        self.er.run()


    def track_vary(self, varydict: dict, varyitemlist=None, mode="row", add_watch_start=False, **kwargs):
        """
        Track a set of particles in combination with a
        vary command.
        """
        assert varyitemlist is not None
        assert len(varyitemlist) == len(varydict)
        assert mode.lower() in ["row", "table"]

        # generate the sdds input file
        sdds = SDDS(self.er.sif, f"{self.er.rootname}.sdds", 0)
        sdds.generate_scan_dataset(varydict)

        self.er.commandfile.clear()
        self.er.add_run_setup()

        if kwargs.get("rf", False):
            self.er.add_rf(volt=kwargs.get("volt", 350000))

        if kwargs.get("rad", False):
            self.er.add_radiation_damping()

        if add_watch_start:
            self.er.add_watch_at_start()
        n_idx = 1 if mode == "row" else len(varydict)
        self.er.commandfile.addCommand(
            "run_control", n_indices=n_idx, n_passes=kwargs.get("n_passes", 2 ** 8)
        )
        if mode == "table":
            for i, it in enumerate(varydict.items()):
                k, v = it
                self.er.add_vary_element_from_file(
                    name=k,
                    item=varyitemlist[i],
                    index_number=i,
                    index_limit=len(v),
                    enumeration_file=f"{self.rootname}.sdds",
                    enumeration_column=k,
                )
        else:
            for i, it in enumerate(varydict.items()):
                k, v = it
                self.er.add_vary_element_from_file(
                    name=k,
                    item=varyitemlist[i],
                    index_number=0,
                    index_limit=len(v),
                    enumeration_file=f"{self.rootname}.sdds",
                    enumeration_column=k,
                )
        self.er.commandfile.addCommand("bunched_beam")
        self.er.commandfile.addCommand(
            "sdds_beam", input=self.er.sdds_beam_file, input_type='"elegant"', reuse_bunch=1
        )
        self.er.commandfile.addCommand("track")
        self.er.run()
        
        
def GenerateNDimCoordinateGrid(N, NPOINTS, pmin=1e-6, pmax=1e-4, man_ranges=None):
    """
    Method to generate an N dimensional coordinate grid for tracking,
    with fixed number of point in each dimension.
    The final shape is printed at creation.

    IMPORTANT:
    Number of grid points scales with N * NPOINTS**N, i.e.
    very large arrays are generated already with
    quite some small numbers for NPOINTS and N.

    Example: NPOINTS = 2, N = 6 -> 6*2*6 = 384 elements

    Parameters:
    -----------
    N: int dimension of the coordinate grid
    NPOINTS: int number of points in each dimension
    pmin: float min coordinate value in each dim
    pmax: float max coordinate value in each dim

    Returns:
    ========
    coordinate_grid : numpy array coordinate grid with particle ID in last column
    """
    rangelist = [np.linspace(pmin, pmax, NPOINTS)] * N
    if man_ranges is not None:
        # print(man_ranges)
        for k, v in man_ranges.items():
            rangelist[int(k)] = v
            # print(rangelist)
    grid = np.meshgrid(*rangelist)
    coordinate_grid = np.array([*grid])
    npart = coordinate_grid.size // N
    coordinate_grid = coordinate_grid.reshape(N, npart).T
    print("Shape: {} - Number of particles: {} ".format(coordinate_grid.shape, npart))
    # add particle id
    coordinate_grid = np.hstack((coordinate_grid, np.array(range(1, npart + 1)).reshape(npart, 1)))
    # print(coordinate_grid)

    return coordinate_grid


def generate_sphere_grid(dim=2, rmin=1e-6, rmax=1, rsteps=3, phisteps=3, **kwargs):
    """Method to generate grid point within n-dim ball, like polar but n-dim.
    Dimension 6 is a special case - as we need it for Elegant tracking. In this case
    the final two dimensions are not polar but fixed for dim 5 and in dim 6 and array
    passed via the kwargs 'deltaGamma'.

    Parameters
    ----------
    dim : int, optional dimension of the ball, by default 2
    rmin : float, optional minimal radius to use, by default 1e-6
    rmax : float, optional maximal radius to use, by default 1
    rsteps : int, optional number of steps in radius grid, by default 3
    phisteps : int, optional number of steps in the angle grid, by default 3
    """
    R = np.linspace(rmin, rmax, rsteps)
    mangle = np.pi

    # only track one kwadrant
    if kwargs.get("half", False):
        mangle = mangle / 2.0

    PHI1 = np.linspace(0, mangle, phisteps)
    PHI2 = np.linspace(0, mangle, phisteps)  # full sphere is 2 pi reduced for tracking to upper half

    # the special case
    if dim != 6:
        matrices = (R,) + tuple((PHI1 for _ in range(dim - 2))) + (PHI2,)
    else:
        # elegant t shift is fixed to zero
        # TODO: fix the fixed t shift
        matrices = (
            (R,)
            + tuple((PHI1 for _ in range(dim - 4)))
            + (PHI2,)
            + (np.array([0.0]), kwargs.get("deltaGamma", np.array([0.0])))
        )

    # create meshgrid to make all combinations
    meshmatrices = np.array(np.meshgrid(*matrices))

    # count the number of particles
    npart = meshmatrices.size // dim

    # reshape
    coord_T = meshmatrices.reshape(dim, npart).T

    #     X = (coord_T[:,0] * np.cos(coord_T[:,1]),)
    X = tuple()

    if dim == 6:
        ndim = 4
    else:
        ndim = dim

    for i in range(1, ndim):
        X += (coord_T[:, 0] * np.prod(np.sin(coord_T[:, 1:i]), axis=1) * np.cos(coord_T[:, i]),)

    X += (coord_T[:, 0] * np.prod(np.sin(coord_T[:, 1:-1]), axis=1) * np.sin(coord_T[:, -1]),)

    if dim != 6:
        sphere_grid = np.vstack(X)
    else:
        sphere_grid = np.vstack(X + (coord_T[:, 4], coord_T[:, 5]))
    print("Shape: {} - Number of paritcles: {} ".format(sphere_grid.T.shape, npart))

    # add particle id
    coordinate_grid = np.hstack((sphere_grid.T, np.array(range(1, npart + 1)).reshape(npart, 1)))
    # print(coordinate_grid)
    return coordinate_grid

        
class particle_generation(ElegantRunToolkit):
    
    def generate_sdds_particle_inputfile(self, grid_type="rectangular", **kwargs):
        """
        Generates an SDDS file containing initial
        particle coordinates on a grid. The grid
        can be defined through the kwargs.

        Parameters:
        ----------
        kwargs      :
        - pmin: min value of grid on each dim
        - pmax: max value of grid on each dim
        - pcentralmev: particle energy (code converts it to beta * gamma )
        - man_ranges: dict containing as key dim num - in order x xp y yp s p and as values an array of values to be used. 
          For p this is autoset to beta gamma based on pcentralmev
        - NPOINTS: number of linear spaced points in each dim for the grid

        Returns:
        --------
        None, writes the data to pre-defined named file.
        """
        assert grid_type in ["rectangular", "spherical"]
        pcentral = kwargs.get("pcentralmev", self.er.energy_gev*1e3)
        print('pcentralmev: ', pcentral)
        # convert to beta * gamma
        pcentral = np.sqrt((pcentral/const.physical_constants["electron mass energy equivalent in MeV"][0])**2  - 1)

        if grid_type == "rectangular":
            npoints_per_dim = kwargs.get("NPOINTS", 2)
            pmin = kwargs.get("pmin", 0)
            pmax = kwargs.get("pmax", 1e-4)
            man_ranges = kwargs.get("man_ranges", {"5": np.array([pcentral])})
            if "5" not in man_ranges.keys() and 5 not in man_ranges.keys():
                man_ranges["5"] = np.array([pcentral])
            # example : man_ranges={'0':np.array([1e-6,1e-5]),'1':[0]})

            # generate coordinate grid, with particle id as last column
            # and save it as plain data table seperated by a whitespace

            gridpoints = GenerateNDimCoordinateGrid(6, npoints_per_dim, pmin=pmin, pmax=pmax, man_ranges=man_ranges)
            particle_df = pd.DataFrame(gridpoints)
            particle_df.to_csv(f"{self.er.rootname}_plain_particles.dat", sep=" ", header=None, index=False)

            # cleanup kwargs
            kwargs.pop("NPOINTS", None)
            kwargs.pop("pmin", None)
            kwargs.pop("pmax", None)
            kwargs.pop("man_ranges", None)
        else:
            rmin = kwargs.get("rmin", 1e-6)
            rmax = kwargs.get("rmax", 1e-1)
            rsteps = kwargs.get("rsteps", 3)
            half = kwargs.get("half", True)
            phisteps = kwargs.get("phisteps", 5)
            deltaGamma = kwargs.get("deltaGamma", np.array([pcentral]))

            particle_df = pd.DataFrame(
                generate_sphere_grid(
                    dim=6,
                    rmin=rmin,
                    rmax=rmax,
                    rsteps=rsteps,
                    phisteps=phisteps,
                    deltaGamma=deltaGamma,
                    half=half,
                )
            )

            particle_df.to_csv(f"{self.er.rootname}_plain_particles.dat", sep=" ", header=None, index=False)
            # clean up kwargs
            kwargs.pop("rmin", None)
            kwargs.pop("rmax", None)
            kwargs.pop("rsteps", None)
            kwargs.pop("half", None)
            kwargs.pop("phisteps", None)
            kwargs.pop("deltaGamma", None)

        kwargs.pop("pcentralmev", None)

        # Create sddscommand object
        sddscommand = SDDSCommand(self.er.sif, rootname=self.er.rootname)

        # update the command parameters
        if self.er.parallel:
            outputmode = "binary"
        else:
            outputmode = "ascii"
        kwargs["outputMode"] = outputmode
        kwargs["file_2"] = (f"{self.er.rootname}_particles_input.txt" if not self.er.parallel else f"{self.er.rootname}_particles_input.bin")

        # load the pre-defined  convert plain data to sdds command
        cmd = sddscommand.get_particles_plain_2_SDDS_command(**kwargs)

        # run the sdds command
        sddscommand.runCommand(cmd)

        sdds_beam_file = kwargs["file_2"]
        return sdds_beam_file

    
        
class dynap(ElegantRunToolkit):
        
    def add_DA_command(self, **kwargs): # TODO: move this to ElegantRun class
        """
        Add DA find aperture command.
        """
        self.er.commandfile.addCommand(
            "find_aperture",
            output="%s.aper",
            mode=kwargs.get("mode", "n-line"),
            verbosity=0,
            xmin=kwargs.get("xmin", -0.1),
            xmax=kwargs.get("xmax", 0.1),
            xpmin=kwargs.get("xpmin", 0.0),
            xpmax=kwargs.get("xpmax", 0.0),
            ymin=kwargs.get("ymin", 0.0),
            ymax=kwargs.get("ymax", 0.1),
            ypmin=kwargs.get("ypmin", 0.0),
            ypmax=kwargs.get("ypmax", 0.0),
            nx=kwargs.get("nx", 21),
            ny=kwargs.get("ny", 11),
            n_lines=kwargs.get("n_lines", 11),
            split_fraction=kwargs.get("split_fraction", 0.5),
            n_splits=kwargs.get("n_splits", 0),
            desired_resolution=kwargs.get("desired_resolution", 0.01),
            offset_by_orbit=kwargs.get("offset_by_orbit", 0),
            full_plane=kwargs.get("full_plane", 1),
        )
        
    def dynap(self, **kwargs):
        """
        Run Elegant's Dynamic Aperture.
        """
        self.er.commandfile.clear()
        self.er.commandfile.addCommand(
            "run_setup",
            lattice=self.lattice,
            use_beamline=self.beamline,
            p_central_mev=self.energy_gev*1e3,
            centroid="%s.cen",
            default_order=kwargs.get("default_order", 2),
            concat_order=kwargs.get("concat_order", 3),
            rootname=self.rootname,
            parameters="%s.params",
            semaphore_file="%s.done",
            magnets="%s.mag",  # for plotting profile
            losses="%s.los",
        )

        self.er.commandfile.addCommand("twiss_output", filename="%s.twi", output_at_each_step=1)
        self.er.commandfile.addCommand("run_control", n_passes=kwargs.pop("n_passes", 2 ** 9))
        self.add_DA_command(**kwargs)

        self.er.run()

    def dynapmom(self):
        """
        Run Elegant's Dynamic Momentum Aperture.
        """
        # TODO
        pass

        
class watchpoints(ElegantRunToolkit):
    
    def add_watch(self, **kwargs):
        """Add watch point."""
        self.er.commandfile.exists('run_setup')
        self.er.commandfile.addCommand(
            "insert_elements",
            name=kwargs.get("name", ""),
            type=kwargs.get("type", ""),
            exclude="",
            s_start=kwargs.get("s_start", -1),
            s_end=kwargs.get("s_end", -1),
            skip=kwargs.get("skip", 1),
            insert_before=kwargs.get("insert_before", 0),
            add_at_end=kwargs.get("add_at_end", 0),
            add_at_start=kwargs.get("add_at_start", 0),
            element_def=kwargs.get(
                "element_def",
                r'"WQ: WATCH, FILENAME=\"%s-%03ld.wq\", mode=\"coordinates\", interval={}"'.format(
                    kwargs.get("interval", 1)
                ),
            ),
        )

    def add_watch_at_start(self):
        """Add watch point at start of lattice."""
        self.add_watch(
            name="W",
            add_at_start=1,
            element_def=r'"W: WATCH, FILENAME=\"%s-%03ld.wq\", mode=\"coordinates\""',
        )


class misc(ElegantRunToolkit):

    def add_vary_element_from_file(self, **kwargs):
        """
        Add single vary element line, loading value from
        dataset file.
        """
        if "enumeration_file" not in kwargs.keys():
            raise RuntimeError("External filename missing.")
        else:
            self.er.commandfile.addCommand(
                "vary_element",
                name=kwargs.get("name", "*"),
                item=kwargs.get("item", "L"),
                index_number=kwargs.get("index_number", 0),
                index_limit=kwargs.get("index_limit", 1),
                enumeration_file=kwargs.get("enumeration_file"),
                enumeration_column=kwargs.get("enumeration_column"),
            )

    def use_standard_nkicks(self):
        self.er.add_alter_elements(name="*", type="CSBEND", item="N_KICKS", value=16)
        self.er.add_alter_elements(name="*", type="KQUAD", item="N_KICKS", value=8)
        self.er.add_alter_elements(name="*", type="KSEXT", item="N_KICKS", value=8)
        
        
    def findtwiss(self, **kwargs):
        """
        Run Twiss and return Twiss parameters
        together with Twiss data.

        Parameters:
        ----------
        kwargs  : dict twiss command options
        """
        # TODO: add matched = 0 case
        matched = kwargs.get("matched", 1)
        initial_optics = kwargs.get("initial_optics", [])
        alternate_element = kwargs.get("alternate_elements", {})
        closed_orbit = kwargs.get("closed_orbit", 1)

        # make sure not residual is there
        self.er.commandfile.clear()

        # add setup command
        self.er.add_run_setup()

        # add twiss calc
        self.er.commandfile.addCommand(
            "twiss_output",
            matched=matched,
            output_at_each_step=0,
            filename="%s.twi",
            radiation_integrals=1,
        )

        # add controls # TODO: not yet working for new version
        self.add_basic_controls()

        # write command file
        self.er.commandfile.write()
        self.er.commandfile.clear()

        # set cmdstr and run
        cmdstr = "{} elegant {}.ele".format(self.er.sif, self.er.rootname)
        with open(os.devnull, "w") as f:
            subp.call(shlex.split(cmdstr), stdout=f)

        # load twiss output
        twifile = SDDS(self.er.sif, "{}.twi".format(self.er.rootname), 0)
        twiparams = twifile.getParameterValues()
        twidata = twifile.getColumnValues()

        twiparams["length"] = np.round(twidata.iloc[-1]["s"], 3)

        return twidata, twiparams

    def find_matrices(self, **kwargs):
        """
        Find element by element matrix and map elements (depending on given order).
        Constant vector and R matrix are returned as numpy arrays, the maps are
        returned as dicts.

        Parameters:
        -----------
        kwargs  :
          - SDDS_output_order : order of maps (max is 3)

        Returns:
        --------
        C       : np.array constant vector
        R       : np.array R matrix
        T_dict  : dict T map Tijk as key
        Q_dict  : dict U map Qijkl as key
        """
        assert kwargs.get("SDDS_output_order", 1) < 4

        self.er.commandfile.clear()
        self.er.add_run_setup()
        self.er.commandfile.addCommand(
            "matrix_output",
            SDDS_output="%s.sdds",
            SDDS_output_order=kwargs.get("SDDS_output_order", 1),
            printout="%s.mat",
            printout_order=kwargs.get("SDDS_output_order", 1),
            full_matrix_only=kwargs.get("full_matrix_only", 0),
            individual_matrices=kwargs.get("individual_matrices", 1),
            output_at_each_step=kwargs.get("output_at_each_step", 1),
        )

        # add controls # TODO: not yet working for new version
        self.add_basic_controls()

        # write command file
        self.er.commandfile.write()
        self.er.commandfile.clear()

        # set cmdstr
        cmdstr = "{} elegant {}.ele".format(self.er.sif, self.er.rootname)
        with open(os.devnull, "w") as f:
            subp.call(shlex.split(cmdstr), stdout=f)

        with open(f"{self.er.rootname}.mat", "r") as f:
            mdata = f.read()

        # get full turn matrix and
        dfmat = pd.read_csv(
            StringIO("\n".join(mdata.split("full", 1)[1].splitlines()[1:])),
            delim_whitespace=True,
            names=[1, 2, 3, 4, 5, 6],
        )
        C = dfmat.loc[dfmat.index == "C:"].values.T
        R = dfmat.loc[dfmat.index.str.contains("R")].values
        T = dfmat.loc[dfmat.index.str.contains("T")]
        Q = dfmat.loc[dfmat.index.str.contains("Q")]

        T_dict = {}
        for _, row in T.iterrows():
            _basekey = row.name[:-1]
            for c in T.columns:
                _key = _basekey + str(c)
                _value = row[c]
                if not pd.isna(_value):
                    T_dict[_key] = _value

        Q_dict = {}
        for _, row in Q.iterrows():
            _basekey = row.name[:-1]
            for c in Q.columns:
                _key = _basekey + str(c)
                _value = row[c]
                if not pd.isna(_value):
                    Q_dict[_key] = _value

        sddsmat = SDDS(self.er.sif, f"{self.er.rootname}.sdds", 0)
        ElementMatrices = sddsmat.getColumnValues()

        return C, R, ElementMatrices, T_dict, Q_dict

    def fma(self, **kwargs):
        """
        Run Elegant frequency map analysis.
        """
        self.er.commandfile.clear()
        self.er.add_run_setup()

        self.er.commandfile.addCommand("run_control", n_passes=kwargs.pop("n_passes", 2 ** 8))
        self.er.add_twiss_output()
        self.er.add_frequency_map(**kwargs)

        self.er.run()

    def table_scan(self, scan_list_of_dicts, mode="row", add_watch_start=True, **kwargs):
        """ """
        assert mode in ["row", "table"]
        n_idx = 1 if mode == "row" else len(scan_list_of_dicts)
        self.er.commandfile.clear()
        self.er.add_run_setup()

        if add_watch_start:
            self.er.add_watch_at_start()
        self.er.commandfile.addCommand(
            "run_control", n_indices=n_idx, n_passes=kwargs.get("n_passes", 2 ** 8)
        )
        for i, l in enumerate(scan_list_of_dicts):
            if mode == "table":
                inx = i
            else:
                inx = 0
            self.er.commandfile.addCommand(
                "vary_element",
                name=l.get("name"),
                item=l.get("item"),
                initial=l.get("initial"),
                final=l.get("final"),
                index_number=inx,
                index_limit=l.get("index_limit"),
            )
        self.er.commandfile.addCommand(
            "sdds_beam", input=self.er.sdds_beam_file, input_type='"elegant"', reuse_bunch=1
        ) # TODO sdds_beam_file handling ...
        self.er.commandfile.addCommand("track")

        # run will write command file and execute it
        self.er.run()
        
        