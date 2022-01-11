from .main import ElegantRunToolkit

from pyelegantsdds.sdds import SDDSCommand

import numpy as np
import pandas as pd
from scipy import constants as const

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
