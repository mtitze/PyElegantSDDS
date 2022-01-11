from .main import ElegantRunToolkit
from .rf import rftools

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

    def dynapmom(self, n_passes, rf=0, **kwargs):
        """
        Run
        Elegant's Dynamic Momentum Aperture.
        """
        self.er.commandfile.clear()
        self.er.add_run_setup(**kwargs)
    
        if rf != 0:
            rfc = rftools(self.er)
            rfc.rf_handling(rf_status=rf, **kwargs)
            
        self.er.commandfile.addCommand("run_control", n_passes=n_passes)
        self.er.add_momentum_aperture(**kwargs)
