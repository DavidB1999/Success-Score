


### *Ball_Control(pos_data, pitch : Pitch, possession_object=None, method: str = 'dfl', distance: int = 1, direction_A1: str = 'ltr', area: str = 'final third', direction_A2 = None)* <br>

**Parameters**: <br>

+ **pos_data** (dict) - Dictionary of floodlight XY-objects providing spatiotemporal data for both teams in both halves; as returned by dfl.read_position_data_xml
+ **pitch** (Pitch)  - A floodlight Pitch object corresponding to the XY data that will be supplied to the model.
+ **method** ({‘dfl’, ‘custom’}, optional) - A string defining the method to be used to determine ball control. With 'dfl' the possession_objects returned by dfl.read_position_data_xml are used to identify ball control events. With 'custom' ball control is defined via the distance of both teams's closest player to the ball.
+ **distance** (int, optional) - The distance in meters used if *method = custom*
+ **direction_A1** ({‘ltr’, ‘rtl’}, optional) - The playing direction for the home team in the first half (left to right or right to left)
+ **area** ({'penalty area', 'final third', '30m', 'half', 'full'}, optional) - Area in which ball control will count as ball control event
+ **directon_A2** ({‘ltr’, ‘rtl’}, optional) - The playing direction for the home team in second half (left to right or right to left). If None, it will be set to the opposite of *direction_A1*

**Returns**: <br>

+ **BC_A1_area, BC_B1_area, BC_A2_area, BC_B2_area** (lists) - List of dummy {0, 1} indicating ball control events for the two teams {Home: A, Away:B} in the two halve



### *Space_Control(pos_data, pitch:Pitch, direction_A1: str = 'ltr', meshtype: str = 'square', area: str = 'final third', resolution: int = 10, direction_A2 = None)*

**Parameters**: <br>

+ **pos_data** (dict) - Dictionary of floodlight XY-objects providing spatiotemporal data for both teams in both halves; as returned by dfl.read_position_data_xml
+ **pitch** (Pitch)  - A floodlight Pitch object corresponding to the XY data that will be supplied to the model.
+ **direction_A1** ({‘ltr’, ‘rtl’}, optional) - The playing direction for the home team in the first half (left to right or right to left)
+ **meshtype** ({‘square’, ‘hexagonal’}, optional) -  A string indicating the type of mesh that will be generated. ‘square’ will generate a grid-like mesh with square cell shapes (default). ‘hexagonal’ will generate a mesh with hexagonal cell shapes where mesh points have equidistant neighbours.
+ **area** ({'penalty area', 'final third', '30m', 'half', 'full'}, optional) - Area in which the controlled spaces will be assessed; for which the meshgrid will be generated
+ **resolution** (int, optional) - The number of mesh grid points used in x-direction. Must be in range [10, 1000] and defaults to 100. The number of mesh grid points in y-direction will be inferred automatically to match the shape of the pitch and produce regular mesh cell shapes.
+ **directon_A2** ({‘ltr’, ‘rtl’}, optional) - The playing direction for the home team in second half (left to right or right to left). If None, it will be set to the opposite of *direction_A1*

**Returns**: <br>

+ **SCR_A1, SCR_B1, SCR_A2, SCR_B2** (TeamProperty) - Space control rates for both teams {Home: A, Away:B} in both halves i.e. One Property object for each team (corresponding to the fitted xy1 and xy2) of shape (n_frames x 1), respectively. Property objects contain the percentage of points controlled by each team on the pitch
+ **dvm1, dvm2** (DiscreteVoronoiModel) -DiscreteVoronoiModel-class-object which can be used for plotting, ...
  


### *Success_Score(position_file, matchinfo, area = 'penalty area', ball_control_method = 'dfl',space_control_resolution = 10, space_control_mesh = 'square', direction_A1 = 'ltr', direction_A2 = None, interval_length = 100, Hz = 1, space_control_rate = 20, ball_control_distance = 1, negatives = False)*

**Parameters**: <br>

+ **filepath_positions** (str or pathlib.Path) – Full path to XML File where the Position data in DFL format is saved.
+ **matchinfo**  (str or pathlib.Path) – Full path to XML File where the Match Information data in DFL format is saved.
+ **area** ({'penalty area', 'final third', '30m', 'half', 'full'}, optional) - Area of interest for the Success-Score - will be passed to *area*-parameter of *Space_Control* and *Ball_Control*
+ **ball_control_method** ({‘dfl’, ‘custom’}, optional) - A string defining the method to be used to determine ball control. Will be passed to *method*-parameter of *Ball_Control*. With 'dfl' the possession_objects returned by dfl.read_position_data_xml are used to identify ball control events. With 'custom' ball control is defined via the distance of both teams' closest player to the ball.
+ **space_control_resolution** (int, optional) - The number of mesh grid points used in x-direction. Will be passed to *resolution*-paramter of *Space_Control*. Must be in range [10, 1000] and defaults to 100. The number of mesh grid points in y-direction will be inferred automatically to match the shape of the pitch and produce regular mesh cell shapes.
+ **space_control_mesh** ({‘square’, ‘hexagonal’}, optional) -  A string indicating the type of mesh that will be generated. Will be passed to *meshtype*-parameter of *Space_Control*. ‘square’ will generate a grid-like mesh with square cell shapes (default). ‘hexagonal’ will generate a mesh with hexagonal cell shapes where mesh points have equidistant neighbours.
+ **direction_A1** ({‘ltr’, ‘rtl’}, optional) - The playing direction for the home team in the first half (left to right or right to left). Will be passed to the *direction_A1*-parameter of both *Space_Control* and *Ball_Control*
+ **directon_A2** ({‘ltr’, ‘rtl’}, optional) - The playing direction for the home team in second half (left to right or right to left). If None, it will be set to the opposite of *direction_A1*.  Will be passed to the *direction_A2*-parameter of both *Space_Control* and *Ball_Control*
+ **interval_length** (int, optional) - The length of the interval used in the calculation of the Success-Score
+ **Hz** (int, optional) - The intended framerate for the positional data. Balance computation time and accuracy. Parameter must be smaller than the framerate of the data and has to be a divisor of the original framerate.
+ **space_control_rate** (int, optional) - The lower limit of space control rates which will be interpreted as space control events.
+ **ball_control_distance** (int, optional) - Relevant only if *ball_control_method = 'custom'*. Defines the distance of the closest players used for the definition of ball control events
+ **negatives** (Boolean, optional) - If False all negative Efficiency values will be neutralized to 0.
  
