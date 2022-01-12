from .main import ElegantRunToolkit
from .sdds import SDDS
import pandas as pd
import numpy as np
import os
from io import StringIO

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
        
        