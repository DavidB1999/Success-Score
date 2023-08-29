from typing import Tuple

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
from scipy.spatial.distance import cdist
from scipy.stats import pointbiserialr as pbr
import copy

from floodlight import XY, Pitch, TeamProperty, PlayerProperty
from floodlight.models.base import BaseModel, requires_fit
from floodlight.io import dfl

# Copyright (c) 2021 floodlight-sports
# credit: https://floodlight.readthedocs.io/en/latest/_modules/floodlight/models/space.html#DiscreteVoronoiModel
# all added code snippets follow "#-#"

class DiscreteVoronoiModel(BaseModel):
    """Calculates discretized versions of the Voronoi tessellation commonly used to
    assess space control.

    Upon instantiation, this model creates a mesh grid that spans the entire pitch with
    a fixed number of mesh points. When calling the
    :func:`~DiscreteVoronoiModel.fit`-method, closest players to the respective mesh
    points are evaluated and their control assigned to players. Thus, cumulative
    controls and controlled areas are calculated on a discretization of the pitch. The
    following calculations can subsequently be queried by calling the corresponding
    methods:

        - Player Space Control --> :func:`~DiscreteVoronoiModel.player_controls`
        - Team Space Control --> :func:`~DiscreteVoronoiModel.team_controls`

    Furthermore, the following plotting methods are available to visualize the model:

        - Plot controlled areas --> :func:`~DiscreteVoronoiModel.plot`
        - Plot mesh grid --> :func:`~DiscreteVoronoiModel.plot_mesh`

    Parameters
    ----------
    pitch: Pitch
        A floodlight Pitch object corresponding to the XY data that will be supplied to
        the model. The mesh created during instantiation will span this pitch.
    mesh: {'square', 'hexagonal'}, optional
        A string indicating the type of mesh that will be generated. 'square' will
        generate a grid-like mesh with square cell shapes (default). 'hexagonal' will
        generate a mesh with hexagonal cell shapes where mesh points have equidistant
        neighbours.
    xpoints: int, optional
        The number of mesh grid points used in x-direction. Must be in range [10, 1000]
        and defaults to 100. The number of messh grid points in y-direction will be
        inferred automatically to match the shape of the pitch and produce regular mesh
        cell shapes.

    Notes
    -----
    The original work by Taki and Hasegawa proposed to use Voronoi tessellations for
    assessing player dominant regions [1]_. This approach has later been simplified by
    using the Euclidean distance when allocating space to players [2]_ , [3]_.
    Instead of computing algebraic Voronoi regions, this model discretizes the problem
    by sampling space control on a finite number of mesh points across the pitch. This
    runs much faster and can be easier to handle. If an appropriate number of mesh
    points is chosen, the resulting error is expected to be negligible given the common
    spatial inaccuracies of tracking data as well as variations in moving players'
    centers of masses.

    References
    ----------
        .. [1] `Taki, T., & Hasegawa, J. (2000). Visualization of dominant region in
            team games and its application to teamwork analysis. Proceedings Computer
            Graphics International 2000, 227–235.
            <https://ieeexplore.ieee.org/document/852338>`_
        .. [2] `Fonseca, S., Milho, J., Travassos, B., & Araújo, D. (2012). Spatial
            dynamics of team sports exposed by Voronoi diagrams. Human Movement
            Science, 31(6), 1652–1659. <https://doi.org/10.1016/j.humov.2012.04.006>`_
        .. [3] `Rein, R., Raabe, D., & Memmert, D. (2017). “Which pass is better?”
            Novel approaches to assess passing effectiveness in elite soccer. Human
            Movement Science, 55, 172–181.
            <https://doi.org/10.1016/j.humov.2017.07.010>`_

    Examples
    --------
    >>> import numpy as np
    >>> from floodlight import XY, Pitch
    >>> from floodlight.models.space import DiscreteVoronoiModel

    >>> # create data and fit model
    >>> xy1 = XY(np.array(((10, 10, 20, 80, 30, 40), (10, 10, np.nan, np.nan, 35, 35))))
    >>> xy2 = XY(np.array(((90, 90, 80, 20, 75, 80), (90, 90, 75, 25, 80, 70))))
    >>> pitch = Pitch.from_template("opta", length=105, width=68)
    >>> dvm = DiscreteVoronoiModel(pitch)
    >>> dvm.fit(xy1, xy2)

    >>> # print player controls [%] for first team
    >>> player_control1, player_control2 = dvm.player_controls()
    >>> print(player_control1.property)
    [[10.63 19.32 21.71]
     [10.35  0.   36.56]]

    >>> # print team controls [%] for first team
    >>> team_control1, team_control2 = dvm.team_controls()
    >>> print(team_control1.property)
    [[51.66]
     [46.91]]
    """

    #-# adding an option (area that allows the space control calculation to be applied to any area of the pitch
    #-# adding a variable that indicates the playing direction (for the home team); important for area-specific Space control calculation
    def __init__(self, pitch: Pitch, mesh: str = "square", xpoints: int = 100, area: str = 'full', home_direction = 'ltr'):
        super().__init__(pitch)

        # input parameter
        self._mesh_type = mesh
        self._xpoints = xpoints
        self._area_type = area
        self._home_direction = home_direction

        # model parameter
        self._meshx_ = None
        self._meshy_ = None
        self._xpolysize_ = None
        self._ypolysize_ = None
        self._T_ = None
        self._N1_ = None
        self._N2_ = None
        self._framerate = None
        self._cell_controls_ = None

        # checks
        valid_mesh_types = ["square", "hexagonal"]
        if mesh not in valid_mesh_types:
            raise ValueError(
                f"Invalid mesh type. Expected one of {valid_mesh_types}, got {mesh}"
            )
        if xpoints < 10 or xpoints > 1000:
            raise ValueError(
                f"Expected xpoints to be in range [10, 1000], got {xpoints}"
            )

        #-# currently supported preprogrammed area choices
        valid_areas = ['full', 'half', '30m', 'penalty area', 'final third']
        if area not in valid_areas:
            raise ValueError(
                f"Invalid area. Expected on of {valid_areas}, got {area}"
            )
            
        #-# check direction
        valid_directions = ['ltr','rtl']
        if home_direction not in valid_directions:
            raise ValueError(
                f'Invalid playing direction. Please supply one of the following options {valid_direction} to home_direction')        
        
        # generate mesh
        self._generate_mesh(mesh, xpoints, area)

    #-# add the option so that the mesh can cover defined areas            
    def _generate_mesh(self, mesh: str = "square", xpoints: int = 100, area: str = 'full', home_direction: str = 'ltr') -> None:
        """Generates a np.meshgrid for a given mesh type."""
        # param
        self._meshx_ = None
        self._meshy_ = None
        self._meshxl_ = None
        self._meshyl_ = None
        self._xpolysize_ = None
        self._ypolysize_ = None
        xmin, xmax = self._pitch.xlim
        xmin_l, xmax_l = self._pitch.xlim
        ymin, ymax = self._pitch.ylim

        #-# adapt the limits of the mesh according to area choice
        #-# these adaptations are based on the dfl pitch limits xlim=(-50.0, 50.0), ylim=(-34.0, 34.0), unit='m'
        #-# and does not (fully) translate to other pitch dimensions and scales
        #-# these limits apply to the full pitch and the right side of the pitch for the smaller areas
        if area == 'full':
                print('Do you really care about the space control on the entire pitch?')
        elif area == 'half':
                xmin = xmax - (xmax - xmin)/2
                xmax_l = xmax_l - (xmax_l - xmin_l)/2
        elif area == '30m':
                xmin = xmax - 30
                xmax_l = xmin_l + 30
        elif area == 'penalty area':
                xmin = xmax - 16.5
                xmax_l = xmin_l + 16.5
                old_ymin = ymin
                ymin = ymax - (ymax - ymin)/2 - 20.16
                ymax = ymax - (ymax - old_ymin)/2 + 20.16
        elif area == 'final third':
                xmin = xmax - (xmax - xmin)*1/3
                xmax_l = xmin_l + (xmax_l - xmin_l) *1/3
        
        
        #-# create mesh for both left and right side of the pitch
        xmaxs = [xmax, xmax_l]
        xmins = [xmin, xmin_l]
        for numerator in [0,1]:
            if mesh == "square":
                # determine square size
                self._xpolysize_ = (xmaxs[numerator] - xmins[numerator]) / xpoints
                self._ypolysize_ = self._xpolysize_

                # derive number of points in y direction
                ypoints = round((ymax - ymin) / self._ypolysize_)
                # re-adjust ypolysize for stretching/rounding in y direction
                self._ypolysize_ = (ymax - ymin) / ypoints

                # get padding
                xpad = self._xpolysize_ * 0.5
                ypad = self._ypolysize_ * 0.5

                # create unilateral and two-dimensional grid points
                x = np.linspace(xmins[numerator] + xpad, xmaxs[numerator] - xpad, xpoints)
                y = np.linspace(ymax - ypad, ymin + ypad, ypoints)
                
                #-# define both the mesh for left and right
                if numerator == 0:
                    self._meshx_, self._meshy_ = np.meshgrid(x, y)
                elif numerator ==1:
                    self._meshxl_, self._meshyl_ = np.meshgrid(x, y)

            elif mesh == "hexagonal":
                # longitudinal spacing of polygons (minus half polygon that's out of bounds)
                xspace = (xmaxs[numerator] - xmins[numerator]) / (xpoints - 0.5)
                # hexagon size (= radius of outer circumcircle)
                self._xpolysize_ = xspace / np.sqrt(3)
                self._ypolysize_ = self._xpolysize_
                # lateral spacing of polygons (by formula)
                yspace = self._xpolysize_ * 1.5
                # longitudinal padding, also offset for odd rows of polygons
                xpad = xspace * 0.5

                # derive number of points in y direction
                ypoints = round((ymax - ymin) / yspace) + 1

                # unilateral and two-dimensional grid points
                x = np.linspace(xmins[numerator], xmaxs[numerator] - xpad, xpoints)
                y = np.linspace(ymax, ymin, ypoints)
                
                #-# define both the mesh for left and right
                if numerator == 0:
                    self._meshx_, self._meshy_ = np.meshgrid(x, y)
                    # add offset for odd rows
                    self._meshx_[1::2, :] += xpad
                elif numerator == 1:
                    self._meshxl_, self._meshyl_ = np.meshgrid(x, y)
                    # add offset for odd rows
                    self._meshxl_[1::2, :] += xpad


            
        

    def _calc_cell_controls(self, xy1: XY, xy2: XY):
        """Calculates xID of closest player to each mesh point at each time point and
        stores results in self._cell_controls"""
        # bin
        T = len(xy1)

        #-# need to calculate space control left and right side
        self._cell_controls_ = np.full(
            # shape is: time x (mesh shape)
            (T, self._meshx_.shape[0], self._meshx_.shape[1]),
            np.nan,
        )
        #-# left
        self._cell_controlsl_ = np.full(
            # shape is: time x (mesh shape)
            (T, self._meshxl_.shape[0], self._meshxl_.shape[1]),
            np.nan,
        )

        # loop
        for t in range(T):
            #-# duplicate code for left side of the pitch
            # stack and reshape player and mesh coordinates to (M x 2) arrays
            player_points = np.hstack((xy1.frame(t), xy2.frame(t))).reshape(-1, 2)
            mesh_points = np.stack((self._meshx_, self._meshy_), axis=2).reshape(-1, 2)
            mesh_pointsl = np.stack((self._meshxl_, self._meshyl_), axis=2).reshape(-1, 2)

            # calculate pairwise distances and determine closest player
            pairwise_distances = cdist(mesh_points, player_points)
            pairwise_distancesl = cdist(mesh_pointsl, player_points)

            closest_player_index = np.nanargmin(pairwise_distances, axis=1)
            closest_player_indexl = np.nanargmin(pairwise_distancesl, axis=1)

            self._cell_controls_[t] = closest_player_index.reshape(self._meshx_.shape)
            self._cell_controlsl_[t] = closest_player_indexl.reshape(self._meshxl_.shape)
      

    def fit(self, xy1: XY, xy2: XY):
        """Fit the model to the given data and calculate control values for mesh points.

        Parameters
        ----------
        xy1: XY
            Player spatiotemporal data of the first team.
        xy2: XY
            Player spatiotemporal data of the second team.
        """
        # derive parameters
        self._N1_ = xy1.N
        self._N2_ = xy2.N
        self._T_ = len(xy1)
        self._framerate = xy1.framerate
        # invoke control calculation
        self._calc_cell_controls(xy1, xy2)
        


    @requires_fit
    def player_controls(self) -> Tuple[PlayerProperty, PlayerProperty]:
        """Returns the percentage of mesh points controlled by each player of the first
        and second team.

        Returns
        -------
        player_controls: Tuple[PlayerProperty, PlayerProperty]
            One Property object for each team (corresponding to the fitted xy1 and xy2)
            of shape (n_frames x n_players), respectively. Property objets contain the
            percentage of points controlled by each player on the pitch.
        """
        #-# duplicate all steps for left (suffix l) and right
        # infer number of mesh cells
        number_of_cells = self._cell_controls_.shape[1] * self._cell_controls_.shape[2]
        number_of_cellsl = self._cell_controlsl_.shape[1] * self._cell_controlsl_.shape[2]


        # xID ranges for both team's players if stacked together
        range1 = range(self._N1_)
        range2 = range(self._N1_, self._N1_ + self._N2_)

        # for each xID count number of cell controls in each mesh through time
        counts1 = [np.sum(self._cell_controls_ == xID, axis=(1, 2)) for xID in range1]
        counts2 = [np.sum(self._cell_controls_ == xID, axis=(1, 2)) for xID in range2]
        
        counts1l = [np.sum(self._cell_controlsl_ == xID, axis=(1, 2)) for xID in range1]
        counts2l = [np.sum(self._cell_controlsl_ == xID, axis=(1, 2)) for xID in range2]

        # transform to arrays and normalize
        counts1 = np.array(counts1).transpose()
        counts2 = np.array(counts2).transpose()
        
        counts1l = np.array(counts1l).transpose()
        counts2l = np.array(counts2l).transpose()

        # transform to percentages
        percentages1 = np.round(100 * counts1 / number_of_cells, 2)
        percentages2 = np.round(100 * counts2 / number_of_cells, 2)
        percentages1l = np.round(100 * counts1l / number_of_cellsl, 2)
        percentages2l = np.round(100 * counts2l / number_of_cellsl, 2)

        # create objects
        property1 = PlayerProperty(
            property=percentages1, name="space control", framerate=self._framerate
        )
        property2 = PlayerProperty(
            property=percentages2, name="space control", framerate=self._framerate
        )

        property1l = PlayerProperty(
            property=percentages1l, name="space control", framerate=self._framerate
        )
        property2l = PlayerProperty(
            property=percentages2l, name="space control", framerate=self._framerate
        )
        return property1, property2, property1l, property2l


    @requires_fit
    def team_controls(self) -> Tuple[TeamProperty, TeamProperty]:
        """Returns the percentage of mesh points controlled by the first and second
        team.

        Returns
        -------
        team_controls: Tuple[TeamProperty, TeamProperty]
            One Property object for each team (corresponding to the fitted xy1 and xy2)
            of shape (n_frames x 1), respectively. Property objets contain the
            percentage of points controlled by each team on the pitch.
        """
        #-# duplicate all steps for left (suffix l) and right

        # infer number of mesh cells
        number_of_cells = self._cell_controls_.shape[1] * self._cell_controls_.shape[2]
        number_of_cellsl = self._cell_controlsl_.shape[1] * self._cell_controlsl_.shape[2]


        # count number of cell controls for a team in each mesh through time
        counts1 = np.sum(self._cell_controls_ < self._N1_, axis=(1, 2))
        counts2 = np.sum(self._cell_controls_ >= self._N1_, axis=(1, 2))
        counts1l = np.sum(self._cell_controlsl_ < self._N1_, axis=(1, 2))
        counts2l = np.sum(self._cell_controlsl_ >= self._N1_, axis=(1, 2))

        # transform to arrays and normalize
        counts1 = np.array(counts1).reshape(-1, 1)
        counts2 = np.array(counts2).reshape(-1, 1)
        counts1l = np.array(counts1l).reshape(-1, 1)
        counts2l = np.array(counts2l).reshape(-1, 1)

        # transform to percentages
        percentages1 = np.round(100 * counts1 / number_of_cells, 2)
        percentages2 = np.round(100 * counts2 / number_of_cells, 2)
        percentages1l = np.round(100 * counts1l / number_of_cellsl, 2)
        percentages2l = np.round(100 * counts2l / number_of_cellsl, 2)

        # create objects
        property1 = TeamProperty(
            property=percentages1, name="space control Team A, right", framerate=self._framerate
        )
        property2 = TeamProperty(
            property=percentages2, name="space control Team B, right", framerate=self._framerate
        )
        property1l = TeamProperty(
            property=percentages1l, name="space control, Team A, left", framerate=self._framerate
        )
        property2l = TeamProperty(
            property=percentages2l, name="space control, Team B, left", framerate=self._framerate
        )

        if self._home_direction == 'ltr':
            return property1, property2l
        elif self._home_direction =='rtl':
            return property1l, property2,

    @requires_fit
    def plot(
        self,
        t: int = 0,
        team_colors: Tuple[str, str] = ("red", "blue"),
        ax: matplotlib.axes = None,
        **kwargs,
    ) -> matplotlib.axes:
        """Plots the fitted mesh grid colored by team controls for a given time point
        on a matplotlib axes.

        Parameters
        ----------
        t: int, optional
            Frame for which controls are plotted. Defaults to 0.
        team_colors: Tuple[str, str], optional
            Tuple of two colors in a format accepted by matplotlib that is used to
            color team specific control areas. Defaults to ('red', 'blue').
        ax: matplotlib.axes, optional
            Axes from matplotlib library to plot on. Defaults to None.
        kwargs:
            Optional keyworded arguments e.g. {'zorder', 'ec', 'alpha'} which can be
            used for the plot functions from matplotlib. The kwargs are only passed to
            all the plot functions of matplotlib. If not given default values are used.

        Returns
        -------
        axes: matplotlib.axes
            Axes from matplotlib library with plot.

        Notes
        -----
        The kwargs are only passed to the plot functions of matplotlib. To customize the
        plots have a look at
        `matplotlib
        <https://matplotlib.org/3.5.0/api/_as_gen/matplotlib.axes.Axes.plot.html>`_.

        Examples
        --------
        Given a DiscreteVoronoiModel that has already been fitted:

        >>> # fitted_dvm_model has square mesh
        >>> ax = pitch.plot(color_scheme="bw")
        >>> fitted_dvm_model.plot(ax=ax)

        .. image:: ../../_img/sample_dvm_plot_square.png

        >>> # fitted_dvm_model has hexagonal mesh
        >>> ax = pitch.plot(color_scheme="bw")
        >>> fitted_dvm_model.plot(ax=ax)

        .. image:: ../../_img/sample_dvm_plot_hex.png
        """
        # get ax
        ax = ax or plt.subplots()[1]

        # get colors and construct team color vector
        team_color1, team_color2 = team_colors
        color_vector = [team_color1] * self._N1_ + [team_color2] * self._N2_

        # call plot by mesh type
        if self._mesh_type == "square":
            ax = self._plot_square(t, color_vector, ax=ax, **kwargs)
        elif self._mesh_type == "hexagonal":
            ax = self._plot_hexagonal(t, color_vector, ax=ax, **kwargs)

        return ax


    def _plot_square(
        self,
        t: int = 0,
        team_colors: Tuple[str, str] = None,
        ax: matplotlib.axes = None,
        **kwargs,
    ) -> matplotlib.axes:
        """Plots square mesh grid controls in given color."""
        # handle kwargs
        ec = kwargs.pop("ec", "grey")
        alpha = kwargs.pop("alpha", 0.3)

        # offset to shift rectangle position from bottom left corner to center
        xoffset = -(self._xpolysize_ * 0.5)
        yoffset = -(self._ypolysize_ * 0.5)
        # loop through mesh points and plot Rectangle patch
        for i, j in np.ndindex(self._meshx_.shape):
            poly = plt.Rectangle(
                (self._meshx_[i, j] + xoffset, self._meshy_[i, j] + yoffset),
                width=self._xpolysize_,
                height=self._ypolysize_,
                fc=team_colors[int(self._cell_controls_[t, i, j])],
                ec=ec,
                alpha=alpha,
                **kwargs,
            )
            ax.add_patch(poly)

        #-# repeat for left side of the pitch    
        for i, j in np.ndindex(self._meshxl_.shape):
            poly = plt.Rectangle(
                (self._meshxl_[i, j] + xoffset, self._meshyl_[i, j] + yoffset),
                width=self._xpolysize_,
                height=self._ypolysize_,
                fc=team_colors[int(self._cell_controlsl_[t, i, j])],
                ec=ec,
                alpha=alpha,
                **kwargs,
            )
            ax.add_patch(poly)
        return ax

    def _plot_hexagonal(
        self,
        t: int = 0,
        team_colors: Tuple[str, str] = None,
        ax: matplotlib.axes = None,
        **kwargs,
    ) -> matplotlib.axes:
        """Plots hexagonal mesh grid controls in given color."""
        # handle kwargs
        ec = kwargs.pop("ec", "grey")
        alpha = kwargs.pop("alpha", 0.3)

        # hexagons are regular polygons with 6 vertices
        n_vertices = 6
        # loop through mesh points and plot RegularPolygon patch
        for (i, j), x in np.ndenumerate(self._meshx_):
            poly = RegularPolygon(
                (x, self._meshy_[i, j]),
                numVertices=n_vertices,
                radius=self._xpolysize_,
                fc=team_colors[int(self._cell_controls_[t, i, j])],
                ec=ec,
                alpha=alpha,
                **kwargs,
            )
            ax.add_patch(poly)

        return ax

    def plot_mesh(self, ax: matplotlib.axes = None) -> matplotlib.axes:
        """Plots the generated mesh on a matplotlib.axes.

        Parameters
        ----------
        ax: matplotlib.axes, optional
            Matplotlib axes on which the mesh points are plotted. If ax is None, a
            default-sized matplotlib.axes object is created.

        Returns
        -------
        axes: matplotlib.axes
            Matplotlib axes on which the mesh points are plotted.

        Examples
        --------
        Given a DiscreteVoronoiModel that has already been fitted:

        >>> ax = pitch.plot(color_scheme="bw")
        >>> fitted_dvm_model.plot_mesh(ax=ax)

        .. image:: ../../_img/sample_dvm_plot_hex_mesh.png
        """
        # get ax
        ax = ax or plt.subplots()[1]
        # plot mesh
        ax.plot(self._meshx_, self._meshy_, "ok", markersize=0.5)

        return ax

    
#-# Everything following was added    
#-# function to return the correct space control rates for both teams in a specific area for both halves:
def Space_Control(pos_data, pitch:Pitch, direction_A1: str = 'ltr',
                  meshtype: str = 'square', area: str = 'final third', resolution: int = 10, direction_A2 = None):

    #-# define data for team and half (based on xy_object structure)
    A1 = pos_data['firstHalf']['Home']
    B1 = pos_data['firstHalf']['Away']
    A2 = pos_data['secondHalf']['Home']
    B2 = pos_data['secondHalf']['Away']

    #-# make sure both halves are associated with the correct playing direction
    if direction_A2 == None:
        if direction_A1 == 'ltr':
            direction_A2 = 'rtl'
        elif direction_A1 == 'rtl':
            direction_A2 = 'ltr'
        else:
            raise ValueError(f'Invalid playing direction in A_direction1. Please supply one of the following options: ["ltr", "rtl"]')

    #-# for the first half:
    dvm1 = DiscreteVoronoiModel(pitch = pitch, mesh = meshtype, xpoints=resolution, area=area, home_direction=direction_A1)
    dvm1.fit(A1, B1)
    SCR_A1, SCR_B1 = dvm1.team_controls()

    #-# for the second half:
    dvm2 = DiscreteVoronoiModel(pitch = pitch, mesh = meshtype, xpoints=resolution, area=area, home_direction=direction_A2)
    dvm2.fit(A2, B2)
    SCR_A2, SCR_B2 = dvm2.team_controls()

    return(SCR_A1, SCR_B1, SCR_A2, SCR_B2, dvm1, dvm2)

#-#
#-# function to determine possession

#-# pos_data requires a dictionary of floodlight xy-objects
#-# pitch requires a floodlight pitch object including axis limits
#-# method: either 'dfl' relying on the possession information supplied in the data or 'custom' based on a distance
#-# distance if 'custom' is chosen as method
#-# direction_A1: either 'ltr' or 'rtl' = playing direction of Team A in first half
#-# area: area of possession for Success-Score
#-# possession object created by floodlight dfl.read_position_data_xml; only required if method = 'dfl'

def Ball_Control(pos_data, pitch : Pitch, possession_object=None, method: str = 'dfl', distance: int = 1, direction_A1: str = 'ltr',
                 area: str = 'final third', direction_A2 = None):
    #-# area specific limits 

    area_codes = {'penalty area': {
                    'left':{
                        'xmin' : pitch.xlim[0],
                        'xmax' : pitch.xlim[0] + 16.5,
                        'ymin' : pitch.ylim[1] - (pitch.ylim[1] - pitch.ylim[0])/2 - 20.16,
                        'ymax' : pitch.ylim[1] - (pitch.ylim[1] - pitch.ylim[0])/2 + 20.16,
                    },
                    'right': {
                        'xmin' : pitch.xlim[1] - 16.5,
                        'xmax' : pitch.xlim[1],
                        'ymin' : pitch.ylim[1] - (pitch.ylim[1] - pitch.ylim[0])/2 - 20.16,
                        'ymax' : pitch.ylim[1] - (pitch.ylim[1] - pitch.ylim[0])/2 + 20.16,
                    }
    },
                  'final third': {
                      'left':{
                          'xmin' : pitch.xlim[0],
                          'xmax' : pitch.xlim[0] + (pitch.xlim[1] - pitch.xlim[0])/3,
                          'ymin' : pitch.ylim[0],
                          'ymax' : pitch.ylim[1],
                        },
                      'right': {
                          'xmin' : pitch.xlim[1] - (pitch.xlim[1] - pitch.xlim[0])/3,
                          'xmax' : pitch.xlim[1],
                          'ymin' : pitch.ylim[0],
                          'ymax' : pitch.ylim[1],
                        }
    }, 
                  '30m': {
                      'left':{
                          'xmin' : pitch.xlim[0],
                          'xmax' : pitch.xlim[0] + 30,
                          'ymin' : pitch.ylim[0],
                          'ymax' : pitch.ylim[1],
                        },
                      'right': {
                          'xmin' : pitch.xlim[1] - 30,
                          'xmax' : pitch.xlim[1],
                          'ymin' : pitch.ylim[0],
                          'ymax' : pitch.ylim[1],
                        }
    },  
                  'half': {
                      'left':{
                          'xmin' : pitch.xlim[0],
                          'xmax' : pitch.xlim[0] + (pitch.xlim[1] - pitch.xlim[0])/2,
                          'ymin' : pitch.ylim[0],
                          'ymax' : pitch.ylim[1],
                        },
                      'right': {
                          'xmin' : pitch.xlim[1] - (pitch.xlim[1] - pitch.xlim[0])/2,
                          'xmax' : pitch.xlim[1],
                          'ymin' : pitch.ylim[0],
                          'ymax' : pitch.ylim[1],
                        }
    }}
    
    BC_A1 = []     #-# storing the ball control events per team and half
    BC_B1 = []
    BC_A2 = []
    BC_B2 = []
    P_b1 = pos_data['firstHalf']['Ball'] #-# ball positions second half
    P_b2 = pos_data['secondHalf']['Ball'] #-# ball
    
    #-# if method is custom we need to identify when teams are in possession
    if method == 'custom':
        P_A1 = pos_data['firstHalf']['Home'] #-# Team A position first half
        P_B1 = pos_data['firstHalf']['Away'] #-# Team B
        P_A2 = pos_data['secondHalf']['Home'] #-# Team A positions second half
        P_B2 = pos_data['secondHalf']['Away'] #-# Team B

    
        #-# for first half
        T = len(P_A1)
        #-# determine closest team for each frame for first half
        #-# storing the closest distances for testing purposes
        A1_dists = []
        B1_dists = []
        A2_dists = []
        B2_dists = []

        for t in range(T):
            #-# reformat coordinates
            A1_points = np.hstack((P_A1.frame(t))).reshape(-1, 2)
            B1_points = np.hstack((P_B1.frame(t))).reshape(-1, 2)
            b1_points = np.hstack((P_b1.frame(t))).reshape(-1, 2)

            #-# calculate each players distance to the ball
            A_distances = cdist(A1_points, b1_points)
            B_distances = cdist(B1_points, b1_points)

            A_closest = np.nanmin(A_distances)
            B_closest = np.nanmin(B_distances)
            A1_dists.append(A_closest)
            B1_dists.append(B_closest)


            if A_closest < distance:
                if B_closest >= distance:
                    BC_A1.append(1)
                    BC_B1.append(0)
                else:
                    BC_A1.append(0)
                    BC_B1.append(0)
            elif B_closest < distance:
                BC_A1.append(0)
                BC_B1.append(1)
            else:
                BC_A1.append(0)
                BC_B1.append(0)        

    
        #-# determine closest team for each frame for second half
        for t in range(T):
            #-# reformat coordinates
            A2_points = np.hstack((P_A2.frame(t))).reshape(-1, 2)
            B2_points = np.hstack((P_B2.frame(t))).reshape(-1, 2)
            b2_points = np.hstack((P_b2.frame(t))).reshape(-1, 2)

            #-# calculate each players distance to the ball
            A_distances = cdist(A2_points, b2_points)
            B_distances = cdist(B2_points, b2_points)

            A_closest = np.nanmin(A_distances)
            B_closest = np.nanmin(B_distances)
            A2_dists.append(A_closest)
            B2_dists.append(B_closest)

            if A_closest < distance:
                if B_closest >= distance:
                    BC_A2.append(1)
                    BC_B2.append(0)
                else:
                    BC_A2.append(0)
                    BC_B2.append(0)
            elif B_closest < distance:
                BC_A2.append(0)
                BC_B2.append(1)
            else:
                BC_A2.append(0)
                BC_B2.append(0)      
                
    elif method == 'dfl':
        #-# use deepcopy to avoid the overwriting in PO 
        PO = copy.deepcopy(possession_object)
        BC_A1 = copy.deepcopy(PO['firstHalf'].code)
        BC_B1 = copy.deepcopy(PO['firstHalf'].code)
        BC_A2 = copy.deepcopy(PO['secondHalf'].code)
        BC_B2 = copy.deepcopy(PO['secondHalf'].code)
        
        #-# convert to 0 and 1 (=possession) for each team and each half
        BC_A1[BC_A1 == 1] = 1
        BC_A1[BC_A1 == 2] = 0

        BC_B1[BC_B1 == 1] = 0
        BC_B1[BC_B1 == 2] = 1

        BC_A2[BC_A2 == 1] = 1
        BC_A2[BC_A2 == 2] = 0

        BC_B2[BC_B2 == 1] = 0
        BC_B2[BC_B2 == 2] = 1
    else:
        raise ValueError(f'Method {method} supplied. Use either "dfl" or "custom"')
    
    #-# at this stage we have BC_A1, BC_B1, BC_A2, BC_B2 from either dfl or custom
    #-# now the location of possession (i.e. the ball) needs to be confirmed to be in the area
    #-# only then we count it as ball control event
    
    #-# ball control booleans for the specific area!
    BC_A1_area = []
    BC_B1_area = []
    BC_A2_area = []
    BC_B2_area = []
    new = [BC_A1_area, BC_B1_area, BC_A2_area, BC_B2_area]

    for v, V in enumerate([BC_A1, BC_B1, BC_A2, BC_B2]):
    
    #-# sort ball position to halves
        if v <= 1:
            ball = P_b1
        elif v == 2 or v ==3:
            ball = P_b2
     
   #-# determine the correct direction based on parameters
        if direction_A1 == 'ltr':                                                  #-# if team A played left to right in the first half
            if direction_A2 == None or direction_A2== 'rtl':                       #-# and switched in the second 
                if v == 0 or v == 3:
                    direction = 'right'
                elif v == 1 or v == 2:
                    direction = 'left'
            elif direction_A2 == 'ltr':                                            #-# and did not switch in the second
                if v == 0 or v == 2:
                    direction = 'right'
                elif v == 1 or v == 3:
                    direction = 'left'                

        if direction_A1 == 'rtl':                                                  #-# if team A played right to left in the first half
            if direction_A2 == None or direction_A2== 'ltr':                       #-# and switched in the second 
                if v == 0 or v == 3:
                    direction = 'left'
                elif v == 1 or v == 2:
                    direction == 'right'
            elif direction_A2 == 'rtl':                                            #-# and did not switch in the second
                if v == 0 or v == 2:
                    direction = 'left'
                elif v == 1 or v == 3:
                    direction == 'right'                


        T = len(V)
        for t in range(T):
            if V[t] == 1:
                x = ball.frame(t)[0] #-# x position of ball
                y = ball.frame(t)[1] #-# y position of ball


                if x <= area_codes[area][direction]['xmax'] and x >= area_codes[area][direction]['xmin'] and y <= area_codes[area][direction]['ymax'] and y < area_codes[area][direction]['ymin']:
                    new[v].append(1)
                else:
                    new[v].append(0)
            else: 
                new[v].append(0)

                
    return BC_A1_area, BC_B1_area, BC_A2_area, BC_B2_area

def Success_Score(position_file, matchinfo, area = 'penalty area', ball_control_method = 'dfl',
                  space_control_resolution = 10, space_control_mesh = 'square', direction_A1 = 'ltr',
                 direction_A2 = None, interval_length = 100, Hz = 1, space_control_rate = 20,
                 ball_control_distance = 1, negatives = False):
    

    #-# read position data and match information file
    position_file = 'PositionData/DFL_04_02_positions_raw_DFL-COM-000001_DFL-MAT-0027AD.xml'
    matchinfo ='MI_Data/DFL_02_01_matchinformation_DFL-COM-000001_DFL-MAT-0027AD-Copy1.xml'

    #-# returns multipe objects; at this point we only care about the positions (xy_objects) and the pitch object (pitch)
    xy_objects, possession_objects, ballstatus_objects, teamsheets, pitch = dfl.read_position_data_xml(filepath_positions=position_file,
                              filepath_mat_info=matchinfo)
    
    #-# determine the scaling factor for the downsizing of the data
    org_framerate = xy_objects['firstHalf']['Home'].framerate
    if Hz != org_framerate:
        if org_framerate % Hz == 0 and xy_objects['firstHalf']['Home'].framerate > Hz:
            factor = int(org_framerate/Hz)
        else:
            raise ValueError(f'Current framerate is {org_framerate} and cannot be converted to {Hz}; please chose a Hz that the framerate ca be divided by')
        print(factor)
        #-# downscale the data to Hz by scalign factor
        xy_object2 = copy.deepcopy(xy_objects)
        xy_object2['firstHalf']['Home'].xy = (xy_object2['firstHalf']['Home'].xy[::factor,:])
        xy_object2['firstHalf']['Away'].xy = (xy_object2['firstHalf']['Away'].xy[::factor,:])
        xy_object2['secondHalf']['Home'].xy = (xy_object2['secondHalf']['Home'].xy[::factor,:])
        xy_object2['secondHalf']['Away'].xy = (xy_object2['secondHalf']['Away'].xy[::factor,:])
    
    else:
        xy_object2 = copy.deepcopy(xy_objects)
        
    
    #-# space control rates
    SCR_A1, SCR_B1, SCR_A2, SCR_B2, dvm1, dvm2 = Space_Control(pos_data=xy_object2, pitch=pitch, area=area,
                                                      resolution= space_control_resolution, meshtype = space_control_mesh,
                                                      direction_A1= direction_A1, direction_A2 = direction_A2)
    
    #-# we only need the values (array)
    SCR_A1 = SCR_A1.property.flatten()
    SCR_B1 = SCR_B1.property.flatten()
    SCR_A2 = SCR_A2.property.flatten()
    SCR_B2 = SCR_B2.property.flatten()
    #-# ball control events
    BC_A1, BC_B1, BC_A2, BC_B2 = Ball_Control(pos_data = xy_object2, pitch = pitch, area='penalty area',
                                                direction_A1= direction_A1, direction_A2 = direction_A2,
                                                method = ball_control_method, possession_object=possession_objects, 
                                               distance = ball_control_distance)
    
    
    #-# if ball control is the dfl vector, we need to downsize those as well
    if Hz != org_framerate and ball_control_method == 'dfl':
        BC_A1 = BC_A1[::factor]
        BC_B1 = BC_B1[::factor]
        BC_A2 = BC_A2[::factor]
        BC_B2 = BC_B2[::factor]

    BC = [BC_A1, BC_B1, BC_A2, BC_B2]
    SCR = [SCR_A1, SCR_B1, SCR_A2, SCR_B2]
    
    #-# store space control rates and ball control events
    SpaceControlRate = {'A1': SCR_A1,
                        'B1': SCR_B1,
                        'A2': SCR_A2,
                        'B2': SCR_B2}
    
    BallControl = {'A1': BC_A1,
                   'B1': BC_B1,
                   'A2': BC_A2,
                   'B2': BC_B2}
    
    
    #-# to store EFFORT per frame for each team in each half; start with a zero for frame 1 (no calculation possible!)
    EFFORT = {'A1': [0],
              'B1': [0],
              'A2': [0],
              'B2': [0]}
    EK = list(EFFORT.keys())
    
    #-# to store EFFICIENCY per frame for each team in each half
    EFFICIENCY = {'A1': [0],
                  'B1': [0],
                  'A2': [0],
                  'B2': [0]}
    
    #-# to store EFFICIENCY per frame for each team in each half
    SUCCESSSCORE = {'A1': [],
                    'B1': [],
                    'A2': [],
                    'B2': []}
    
    #-# for each team in each half
    for th,TH in enumerate(BC): 
        
        L = len(TH)                       #-# number of frames
        IL = interval_length * Hz         #-# length of interval in number of frames
        
        #-# EFFORT
        
        #-# for each frame in that half for that team (t = Point of Success-Score)
        for t in range(1, L):
            E = 0             #-# Effort or-sum
            if t < IL:        #-# timepoint of calculation less than 1 interval length away from zero?
                T0 = 0        #-# start calculation at 0 (i.e. shorter interval)
                IL2 = t       #-# means our interval is only as long as the value of t
            else:
                T0 = t - IL   #-# start at t0 of the interval
                IL2 = IL      #-# interval is as long as it was supposed to be

            #-# for each frame in the interval from T0 to t 
            E = 0
            for T in range(T0, t):
                if BC[th][T] == 1 or SCR[th][T] >= space_control_rate:       #-# counts to effort if either ball control or space control rate > value
                    E = E+1
            
            EFFORT[EK[th]].append(E/IL2)                  # append Effort for each Interval / t | Effort = or-sum / interval length
        
        #-# EFFICIENCY
        
        #-# for each frame in that half for that team (t = Point of Success-Score)
        for t in range(1, L):
            
            if t < IL:                               #-# timepoint of calculation less than 1 interval length away from zero?
                T0 = 0                               #-# start calculation at 0 (i.e. shorter interval)
            else:
                T0 = t - IL                          #-# start at t0 of the interval
            BC_int = BC[th][T0:T]                    #-# Ball control for each frame in that interval
            #print(BC_int)
            SCR_int = SCR[th][T0:T]        #-# Space control rate for each frame in that interval
            #print(SCR_int)
            r, p = pbr(BC_int, SCR_int)
            EFFICIENCY[EK[th]].append(r)  

        #-# SUCCESS-SCORE
        
        #-# it is a fair assumption that in almost all cases an interval full of identical ball control events (either 0 or 1)
        #-# are zeros (100 seconds ball control in area less likely than 100 seconds no ball control in area) 
        #-# so we make it zero if it is not calculatable
        EFFICIENCY[EK[th]] = np.nan_to_num(EFFICIENCY[EK[th]], nan=0)
        #-# if we work under the assumption that negative correlation is no worse than no correlation at all we neutralize all negative correlations
        if negatives == False:
            EFFICIENCY[EK[th]][EFFICIENCY[EK[th]]<0] = 0

        for t in range(0, L):
            SUCCESSSCORE[EK[th]].append(EFFICIENCY[EK[th]][t] * EFFORT[EK[th]][t])
            
    return(SUCCESSSCORE, EFFORT, EFFICIENCY, SpaceControlRate, BallControl)
