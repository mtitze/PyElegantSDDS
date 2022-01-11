from .main import ElegantRunToolkit

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
