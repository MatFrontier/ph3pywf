# coding: utf-8

__author__ = "Kerui Lai"
__email__ = "kerui.lai@mail.mcgill.ca"

from atomate.vasp.database import VaspCalcDb
from pymatgen.core import Structure
import numpy as np
import matplotlib.pyplot as plt
from pymatgen.phonon.bandstructure import PhononBandStructureSymmLine
from pymatgen.phonon.dos import CompletePhononDos
from pymatgen.phonon.plotter import PhononBSPlotter, PhononDosPlotter
import csv


class Ph3py_Result:
    def __init__(
        self,
        task_label,
        path_to_db_json,
    ):
        self.has_fitted = False

        # connect to DB
        atomate_db = VaspCalcDb.from_db_file(path_to_db_json)

        # use the Ph3py collection
        ph3py_coll = atomate_db.db["ph3py_tasks"]

        # get our task
        self.ph3py_entry = ph3py_coll.find_one({"task_label": task_label})
        if bool(self.ph3py_entry):
            print(f"Unable to find doc from DB using {task_label = }")
            raise

        self.formula_pretty = self.ph3py_entry["formula_pretty"]

        # set up the pymatgen structure
        self.struct_unitcell = Structure.from_dict(self.ph3py_entry["structure"])

        # extract thermal conductivity calculation result (initial)
        self.T_calc = np.array(self.ph3py_entry["temperature"])
        self.kappa_calc = np.array(self.ph3py_entry["kappa"])
        self.kappa_calc_iso = np.mean(self.kappa_calc[:, :3], axis=1)

        # extract thermal conductivity calculation result (fitted) if available
        if "temperature_fitted" in self.ph3py_entry:
            self.has_fitted = True
            self.T_calc_fitted = np.array(self.ph3py_entry["temperature_fitted"])
            self.kappa_calc_fitted = np.array(self.ph3py_entry["kappa_fitted"])
            self.kappa_calc_iso_fitted = np.mean(self.kappa_calc_fitted[:, :3], axis=1)

        # extract band structure
        self.bs = PhononBandStructureSymmLine.from_dict(
            self.ph3py_entry["band_structure"]
        )

        # extract dos
        self.dos = CompletePhononDos.from_dict(self.ph3py_entry["dos"])

    def plot_thermal_conductivity(
        self,
        ref_filenames=[],
        ref_labels=[],
        plot_dircs=False,
        plot_fitted=True,
        save_file=True,
    ):
        # plot params
        plt.rcParams["font.family"] = "Times New Roman"
        plt.rcParams["font.size"] = 12

        # get reference results
        T_exp = []
        kappa_exp = []
        for filename in ref_filenames:
            T_exp.append([])
            kappa_exp.append([])
            with open(filename, newline="") as f:
                reader = csv.reader(f)
                for row in reader:
                    T_exp[-1].append(float(row[0]))
                    kappa_exp[-1].append(float(row[1]))

        # plot directional calculated kappa (anisotropic)
        if plot_dircs:
            for dirc in range(3):
                plt.plot(
                    self.T_calc,
                    self.kappa_calc[:, dirc],
                    "b:",
                    label=f"Calculated dir{dirc} (initial)",
                )

        # plot calculated kappa (initial)
        plt.plot(
            self.T_calc,
            self.kappa_calc_iso,
            "b--",
            label="Calculated isotropic thermal conductivity (initial)",
        )

        # plot fitted
        if plot_fitted and self.has_fitted:
            plt.plot(
                self.T_calc_fitted,
                self.kappa_calc_iso_fitted,
                "b-",
                label="Calculated isotropic thermal conductivity (fitted) (monoclinic)",
            )

        # plot experiment results
        for i in range(len(self.T_exp)):
            plt.scatter(self.T_exp[i], self.kappa_exp[i], label=ref_labels[i])

        # figure settings
        plt.title("Thermal conductivity of {}".format(self.formula_pretty))
        # plt.xscale("log") # LOG SCALE
        # plt.yscale("log") # LOG SCALE
        plt.xlabel("T (K)")
        plt.ylabel("Thermal conductivity (W/m-K)")
        # plt.xlim(left=200)
        plt.ylim(bottom=0, top=None)
        plt.legend()

        # save the figure
        if save_file:
            plt.savefig(
                "{}-heat-conductivity-{}.png".format(
                    self.formula_pretty, self.task_label
                )
            )

        # show plot
        plt.show()

    def plot_bs(
        self,
        save_file=True,
    ):
        plotter = PhononBSPlotter(self.bs)
        plotter.show()
        if save_file:
            plotter.save_plot(
                "{}-phonon-dispersion-{}.png".format(
                    self.formula_pretty, self.task_label
                ),
                "png",
            )

    def plot_dos(
        self,
        save_file=True,
    ):
        plotter = PhononDosPlotter()
        plotter.add_dos("Total DOS", self.dos)
        plotter.add_dos_dict(self.dos.get_element_dos())
        plotter.show(units="cm-1")
        if save_file:
            plotter.save_plot(
                "{}-phonon-DOS-{}.png".format(self.formula_pretty, self.task_label),
                "png",
                units="cm-1",
            )
