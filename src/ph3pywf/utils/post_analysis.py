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
import itertools


class Ph3py_Result:
    def __init__(
        self,
        task_label,
        path_to_db_json,
    ):
        self.has_fitted = False
        self.task_label = task_label

        # connect to DB
        atomate_db = VaspCalcDb.from_db_file(path_to_db_json)

        # use the Ph3py collection
        ph3py_coll = atomate_db.db["ph3py_tasks"]

        # get our task
        self.ph3py_entry = ph3py_coll.find_one({"task_label": task_label})
        if not bool(self.ph3py_entry):
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
        plot_initial=True,
        plot_dircs=False,
        plot_fitted=True,
        save_file=True,
        fig_size=(8, 6),
        title=None,
        ymax=None,
    ):
        # plot params
        plt.rcParams["font.family"] = "Times New Roman"
        plt.rcParams["font.size"] = 12
        colors = itertools.cycle(
            ["r", "c", "b", "g", "y", "m"]
        )  # color sequence for exp
        markers = itertools.cycle(
            ["o", "s", "v", "^", "D", "P"]
        )  # marker sequence for exp

        # initialize figure
        plt.figure(figsize=fig_size, facecolor="#FFFFFF")

        # get reference results
        T_ref = []
        kappa_ref = []
        for filename in ref_filenames:
            T_ref.append([])
            kappa_ref.append([])
            with open(filename, newline="") as f:
                reader = csv.reader(f)
                for row in reader:
                    T_ref[-1].append(float(row[0]))
                    kappa_ref[-1].append(float(row[1]))

        # plot directional calculated kappa (anisotropic)
        if plot_dircs:
            for dirc in range(3):
                plt.plot(
                    self.T_calc,
                    self.kappa_calc[:, dirc],
                    "k:",
                    label=f"Calculated dir{dirc} (initial)",
                )

        # plot calculated kappa (initial)
        if plot_initial:
            plt.plot(
                self.T_calc,
                self.kappa_calc_iso,
                "k-",
                # label="Calculated isotropic thermal conductivity (initial)",
                label="This work",
            )

        # plot calculated kappa (fitted)
        if plot_fitted and self.has_fitted:
            plt.plot(
                self.T_calc_fitted,
                self.kappa_calc_iso_fitted,
                "k-",
                # label="Calculated isotropic thermal conductivity (fitted)",
                label="This work",
            )

        # plot experiment results
        for i in range(len(T_ref)):
            plt.scatter(
                T_ref[i],
                kappa_ref[i],
                c=next(colors),
                marker=next(markers),
                label=ref_labels[i],
            )

        # figure settings
        if title:
            plt.title(title)
        # plt.title("Thermal conductivity of {}".format(self.formula_pretty))
        # plt.xscale("log") # LOG SCALE
        # plt.yscale("log") # LOG SCALE
        plt.xlabel("Temperature (K)")
        plt.ylabel(r"$\kappa$ (W.$m^{-1}.K^{-1}$)")
        plt.xlim(left=200)
        plt.ylim(bottom=0, top=ymax)
        plt.legend()
        plt.grid(visible=True, linestyle=":")

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
        plotter.show(units="THz")
        if save_file:
            plotter.save_plot(
                "{}-phonon-bs-{}.png".format(self.formula_pretty, self.task_label),
                "png",
                units="THz",
            )

    def plot_dos(
        self,
        save_file=True,
    ):
        plotter = PhononDosPlotter()
        plotter.add_dos("Total DOS", self.dos)
        plotter.add_dos_dict(self.dos.get_element_dos())
        plotter.show(units="THz")
        if save_file:
            plotter.save_plot(
                "{}-phonon-DOS-{}.png".format(self.formula_pretty, self.task_label),
                "png",
                units="THz",
            )
