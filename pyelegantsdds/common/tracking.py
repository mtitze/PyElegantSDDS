from .main import ElegantRunToolkit
from .observe import watchpoints
from .rf import rftools
from .radiation import radiation

import numpy as np

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
        
        if rf != 0:
            rfc = rftools(self.er)
            rfc.rf_handling(rf_status=rf, **kwargs)
            
        if rad:
            radc = radiation(self.er)
            radc.add_radiation_damping(**kwargs)
            
        self.er.commandfile.addCommand("run_control", n_passes=n_passes)

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
