"""

Color is tuple of (R, G, B) where RGB is from 0.0 to 1.0.



Functions:
scatter
line
bar
histogram

calc_histogram

plot_heatmap
plot_heatmap_cmd
find_tall_heatmap_size
find_wide_heatmap_size

place_ticks
position_labels

"""
# _make_graph
# _set_default_color
# _set_default_point_size
# _set_default_shape
# _set_default_label
# _set_default_error_bar
# _set_default_axes    lim, tick, tick_label
# _set_default_lim
# _set_default_tick
# _set_default_tick_label
# _set_bar_width
# 
# _choose_tick_delta
#
# _overlaps
# _overlaps_1d
# _overlaps_3d
# _score_label
# _calc_label_coord
# _position_one_label

class Graph:
    # Assume a drawing library using the 4th quadrant in the Cartesian
    # coordinate system. (0, 0, 0) is at the upper left corner of the
    # screen, extending positive to the right, down, and toward the
    # viewer.
    #
    # Color is tuple of (R, G, B) where RGB is from 0.0 to 1.0.
    #
    # Members:
    # TOTAL_WIDTH   Number of pixels for the whole canvas.
    # TOTAL_HEIGHT
    # TOTAL_DEPTH
    #
    # GRAPH_X       The pixel coordinate where the drawable area starts.
    # GRAPH_Y
    # GRAPH_Z
    # GRAPH_WIDTH
    # GRAPH_HEIGHT
    # GRAPH_DEPTH
    #
    # UNIT          A standard unit for setting relative sizes of things.
    # 
    # _plotter
    # _image
    # _virt_xlim
    # _virt_ylim
    # _virt_zlim
    # _zstack
    #
    # Methods:
    # draw_axes
    # draw_tick_marks
    # draw_labels
    # 
    # draw_points
    # draw_line
    # draw_bars
    # 
    # write
    # 
    # _virtx2pix
    # _virty2pix
    def __init__(
        self, plotter, pixel_width, pixel_height,
        virtual_xlim, virtual_ylim, virtual_zlim):
        # Keeps track of two sets of coordinates:
        # pixels   How things are actually drawn.
        # virtual  Virtual units on the graph.
        
        # Constants.
        self.TOTAL_WIDTH = pixel_width
        self.TOTAL_HEIGHT = pixel_height
        self.TOTAL_DEPTH = max(pixel_width, pixel_height)

        # Calculate the borders around the plot.
        #BORDER = 0.25*min(self.TOTAL_WIDTH, self.TOTAL_HEIGHT)
        XBORDER = 0.25 * self.TOTAL_WIDTH
        YBORDER = 0.25 * self.TOTAL_HEIGHT
        X1_BORDER = XBORDER * 0.65   # Make a little larger for the labels.
        X2_BORDER = XBORDER - X1_BORDER
        Y1_BORDER = YBORDER * 0.35
        Y2_BORDER = YBORDER - Y1_BORDER

        # Upper left hand coordinate of the graph.
        self.GRAPH_X = X1_BORDER
        self.GRAPH_Y = Y1_BORDER
        self.GRAPH_Z = 0
        self.GRAPH_WIDTH = self.TOTAL_WIDTH-(X1_BORDER+X2_BORDER)
        self.GRAPH_HEIGHT = self.TOTAL_HEIGHT-(Y1_BORDER+Y2_BORDER)
        self.GRAPH_DEPTH = self.TOTAL_DEPTH
        
        self.UNIT = max(self.GRAPH_WIDTH, self.GRAPH_HEIGHT)/100.0

        is_3d = virtual_zlim is not None
        self._plotter = plotter
        self._image = plotter.image(
            self.TOTAL_WIDTH, self.TOTAL_HEIGHT, self.TOTAL_DEPTH, is_3d)
        self._virt_xlim = virtual_xlim
        self._virt_ylim = virtual_ylim
        self._virt_zlim = virtual_zlim

        self._zstack = []  # list of (object, z, depth)

    def draw_axes(self, draw_x, draw_y, draw_z, draw_box):
        import graphconst as gc
        
        # Axes are half cylinders lying on the background.
        AXIS_RADIUS = 0.2*self.UNIT
        AXIS_COLOR = 0.4, 0.4, 0.4

        x_min, y_min = self.GRAPH_X, self.GRAPH_Y+self.GRAPH_HEIGHT
        x_axis_at, y_axis_at = y_min, x_min

        if draw_box:
            draw_x = draw_y = draw_z = True
        is_3d = bool(self._virt_zlim)


        # Draw the X axis (the bottom line of the X-Y plane).
        coord = [x_min, x_axis_at, self.GRAPH_Z]
        extent = [self.GRAPH_WIDTH, 0, 0]
        if draw_x:
            self._plotter.cylinder(
                self._image, coord, extent, AXIS_RADIUS, AXIS_COLOR,
                finish=gc.METALLIC, blob=True)
        if draw_box and is_3d:
            coord[2] = self.GRAPH_Z + self.GRAPH_DEPTH
            self._plotter.cylinder(
                self._image, coord, extent, AXIS_RADIUS, AXIS_COLOR,
                finish=gc.METALLIC, blob=True)
            coord[2] = self.GRAPH_Z
        # Draw the top line of the X-Y plane.
        coord[1] = self.GRAPH_Y
        if draw_box:
            self._plotter.cylinder(
                self._image, coord, extent, AXIS_RADIUS, AXIS_COLOR,
                finish=gc.METALLIC, blob=True)
        if draw_box and is_3d:
            coord[2] = self.GRAPH_Z + self.GRAPH_DEPTH
            self._plotter.cylinder(
                self._image, coord, extent, AXIS_RADIUS, AXIS_COLOR,
                finish=gc.METALLIC, blob=True)

        # Draw the Y axis (the left line of the X-Y plane).
        coord = [y_axis_at, y_min, self.GRAPH_Z]
        extent = [0, -self.GRAPH_HEIGHT, 0]
        if draw_y:
            self._plotter.cylinder(
                self._image, coord, extent, AXIS_RADIUS, AXIS_COLOR,
                finish=gc.METALLIC, blob=True)
        if draw_box and is_3d:
            coord[2] = self.GRAPH_Z + self.GRAPH_DEPTH
            self._plotter.cylinder(
                self._image, coord, extent, AXIS_RADIUS, AXIS_COLOR,
                finish=gc.METALLIC, blob=True)
            coord[2] = self.GRAPH_Z
        # Draw the right line of the X-Y plane.
        coord[0] = self.GRAPH_X + self.GRAPH_WIDTH
        if draw_box:
            self._plotter.cylinder(
                self._image, coord, extent, AXIS_RADIUS, AXIS_COLOR,
                finish=gc.METALLIC, blob=True)
        if draw_box and is_3d:
            coord[2] = self.GRAPH_Z + self.GRAPH_DEPTH
            self._plotter.cylinder(
                self._image, coord, extent, AXIS_RADIUS, AXIS_COLOR,
                finish=gc.METALLIC, blob=True)

        # Draw the Z axis (the left line of the X-Z plane).
        # Z-axis.
        coord = [y_axis_at, x_axis_at, self.GRAPH_Z]
        extent = [0, 0, self.GRAPH_DEPTH]
        if draw_z:
            self._plotter.cylinder(
                self._image, coord, extent, AXIS_RADIUS, AXIS_COLOR,
                finish=gc.METALLIC, blob=True)
        if draw_box and is_3d:
            coord[1] = self.GRAPH_Y
            self._plotter.cylinder(
                self._image, coord, extent, AXIS_RADIUS, AXIS_COLOR,
                finish=gc.METALLIC, blob=True)
            coord[1] = self.GRAPH_Y+self.GRAPH_HEIGHT
        # Draw the right line of the X-Z plane.
        coord[0] = self.GRAPH_X + self.GRAPH_WIDTH
        if draw_box and is_3d:
            self._plotter.cylinder(
                self._image, coord, extent, AXIS_RADIUS, AXIS_COLOR,
                finish=gc.METALLIC, blob=True)
            coord[1] = self.GRAPH_Y
            self._plotter.cylinder(
                self._image, coord, extent, AXIS_RADIUS, AXIS_COLOR,
                finish=gc.METALLIC, blob=True)

    def draw_tick_marks(
        self, xtick, ytick, ztick, vertical_xtick_label,
        xtick_label, ytick_label, ztick_label,
        xtick_at, ytick_at, ztick_at, draw_grid, 
        tick_size=1.0, label_size=1.0):
        # xtick is a list of virtual units where the tickmarks should
        # be on the X-axis.
        # wrong_ztick_label will draw the ztick label on the opposite
        # side of the graph.
        import graphconst as gc

        # Tick marks are half cylinders lying on the background.
        TICK_COLOR = 0.4, 0.4, 0.4
        #TICK_LENGTH = 1.6 * self.UNIT
        TICK_LENGTH = 0.8 * self.UNIT
        TICK_RADIUS = 0.2*self.UNIT*tick_size

        GRID_COLOR = 0.8, 0.8, 0.8
        GRID_RADIUS = TICK_RADIUS * 0.50

        LABEL_COLOR = 0.4, 0.4, 0.4
        LABEL_FONTSIZE = 1.5*self.UNIT*label_size
        LABEL_DEPTH = 0.1*self.UNIT
        #LABEL_DEPTH = 0.02*self.UNIT
        LABEL_SPACE = 0.8*self.UNIT   # between tick mark and label

        x_min, y_min = self.GRAPH_X, self.GRAPH_Y+self.GRAPH_HEIGHT
        x_axis_at, y_axis_at = y_min, x_min

        # Convert xtick and ytick to pixels.
        xtick_pixel = [self._virtx2pix(x) for x in xtick]
        ytick_pixel = [self._virty2pix(y) for y in ytick]
        ztick_pixel = []
        if ztick:
            ztick_pixel = [self._virtz2pix(z) for z in ztick]

        for i, x_pix in enumerate(xtick_pixel):
            #coord = [x_pix, x_axis_at-TICK_LENGTH/2.0, 0]
            coord = [x_pix, x_axis_at, 0]
            extent = 0, TICK_LENGTH, 0
            if not xtick_at&gc.BACK:
                coord[2] = self.GRAPH_Z + self.GRAPH_DEPTH - LABEL_DEPTH/2.0
            self._plotter.cylinder(
                self._image, coord, extent, TICK_RADIUS, TICK_COLOR,
                finish=gc.METALLIC, blob=True)
            # Round out the ticks so they look smoother.
            #self._plotter.sphere(
            #    self._image, coord, TICK_RADIUS, TICK_COLOR,
            #    finish=gc.METALLIC)
            #coord = coord[0], coord[1]+TICK_LENGTH, coord[2]
            #self._plotter.sphere(
            #    self._image, coord, TICK_RADIUS, TICK_COLOR,
            #    finish=gc.METALLIC)

            #if draw_grid and ztick:
            if draw_grid:
                coord = x_pix, x_axis_at, 0
                # Draw the grid lines on the bottom of the plot.
                extent = 0, 0, self.GRAPH_DEPTH
                self._plotter.cylinder(
                    self._image, coord, extent, GRID_RADIUS, GRID_COLOR,
                    finish=gc.METALLIC, blob=True)
                # Draw the grid on the back pane of the plot.
                #extent = 0, -self.GRAPH_HEIGHT, 0
                #self._plotter.cylinder(
                #    self._image, coord, extent, TICK_RADIUS, TICK_COLOR,
                #    finish=self._plotter.gc.METALLIC)

            # Draw the tick labels.
            if not xtick_label:
                continue
            if i == 0 and x_pix == y_axis_at and xtick_at&gc.TOP:
                # don't overlap with axis
                continue
            label = xtick_label[i]
            if label:
                # text coord left, top, back.  Extends to the right,
                # down, to user.
                coord = [x_pix, x_axis_at+TICK_LENGTH/2.0+LABEL_SPACE, 0]
                wrong_y = False
                if xtick_at&gc.TOP:
                    coord[1] = x_axis_at-TICK_LENGTH/2.0-LABEL_SPACE
                    wrong_y = True
                if not xtick_at&gc.BACK:
                    coord[2] = self.GRAPH_Z+self.GRAPH_DEPTH-LABEL_DEPTH/2.0
                self._plotter.text(
                    self._image, label, coord, LABEL_DEPTH, LABEL_FONTSIZE,
                    LABEL_COLOR, center_x=True, wrong_y=wrong_y,
                    vertical=vertical_xtick_label)
            
        for i, y_pix in enumerate(ytick_pixel):
            #coord = [y_axis_at-TICK_LENGTH/2.0, y_pix, 0]
            coord = [y_axis_at-TICK_LENGTH, y_pix, 0]
            extent = TICK_LENGTH, 0, 0
            if not ytick_at&gc.BACK:
                coord[2] = self.GRAPH_Z + self.GRAPH_DEPTH - LABEL_DEPTH/2.0
            self._plotter.cylinder(
                self._image, coord, extent, TICK_RADIUS, TICK_COLOR,
                finish=gc.METALLIC, blob=True)
            # Round out the ticks so they look smoother.
            #self._plotter.sphere(
            #    self._image, coord, TICK_RADIUS, TICK_COLOR,
            #    finish=gc.METALLIC)
            ##coord = coord[0]+TICK_LENGTH, coord[1], coord[2]
            ##self._plotter.sphere(
            ##    self._image, coord, TICK_RADIUS, TICK_COLOR,
            ##    finish=gc.METALLIC)

            if draw_grid:
                # Draw the grid on the X-Y plane.
                coord = y_axis_at, y_pix, 0
                extent = self.GRAPH_WIDTH, 0, 0
                self._plotter.cylinder(
                    self._image, coord, extent, GRID_RADIUS, GRID_COLOR,
                    finish=gc.METALLIC, blob=True)
            if draw_grid and ztick:
                # Draw the grid on the Y-Z plane.
                coord = y_axis_at, y_pix, 0
                extent = 0, 0, self.GRAPH_DEPTH
                self._plotter.cylinder(
                    self._image, coord, extent, GRID_RADIUS, GRID_COLOR,
                    finish=gc.METALLIC, blob=True)
            
            if not ytick_label:
                continue
            # Don't overlap with the Z label.
            if i == 0 and ztick_label and ztick_at&gc.LEFT:
                continue
            label = ytick_label[i]
            if label:
                coord = [y_axis_at-TICK_LENGTH/2.0-LABEL_SPACE, y_pix, 0]
                if not ytick_at&gc.BACK:
                    coord[2] = self.GRAPH_Z+self.GRAPH_DEPTH-LABEL_DEPTH/2.0
                self._plotter.text(
                    self._image, label, coord, LABEL_DEPTH, LABEL_FONTSIZE,
                    LABEL_COLOR, wrong_x=True, center_y=True)

        for i, z_pix in enumerate(ztick_pixel):
            #coord = [y_axis_at-TICK_LENGTH/2.0, x_axis_at, z_pix]
            #extent = TICK_LENGTH, 0, 0
            coord = [y_axis_at-TICK_LENGTH, x_axis_at, z_pix]
            extent = TICK_LENGTH, 0, 0
            if not ztick_at&gc.LEFT:
                coord[0] += self.GRAPH_WIDTH
            self._plotter.cylinder(
                self._image, coord, extent, TICK_RADIUS, TICK_COLOR,
                finish=gc.METALLIC, blob=True)
            # Round out the ticks so they look smoother.
            #self._plotter.sphere(
            #    self._image, coord, TICK_RADIUS, TICK_COLOR,
            #    finish=gc.METALLIC)
            ##coord = coord[0]+TICK_LENGTH, coord[1], coord[2]
            ##self._plotter.sphere(
            ##    self._image, coord, TICK_RADIUS, TICK_COLOR,
            ##    finish=gc.METALLIC)

            if draw_grid:
                # Draw the grid on the X-Z plane.
                coord = y_axis_at, x_axis_at, z_pix
                extent = self.GRAPH_WIDTH, 0, 0
                self._plotter.cylinder(
                    self._image, coord, extent, GRID_RADIUS, GRID_COLOR,
                    finish=gc.METALLIC, blob=True)

            if not ztick_label:
                continue
            label = ztick_label[i]
            if label:
                # text coord left, top, back.  Extends to the right,
                # down, to user.
                coord = [
                    y_axis_at-TICK_LENGTH/2.0-LABEL_SPACE, x_axis_at, z_pix]
                wrong_x = True
                if not ztick_at&gc.LEFT:
                    coord[0] = y_axis_at + self.GRAPH_WIDTH + LABEL_SPACE
                    wrong_x = False
                self._plotter.text(
                    self._image, label, coord, LABEL_DEPTH, LABEL_FONTSIZE,
                    LABEL_COLOR, wrong_x=wrong_x, center_y=True, center_z=True)

    def draw_title(self, title, title_size=1.0):
        # Title sits on the background.
        TITLE_DEPTH = 0.8*self.UNIT
        TITLE_COLOR = 0.2, 0.2, 0.2
        TITLE_FONTSIZE = 4.0*self.UNIT*title_size
        TITLE_SPACE = 0.1*self.UNIT   # space between graph and title

        if not title:
            return
        
        GRAPH_X_MID = self.GRAPH_X+self.GRAPH_WIDTH/2.0
        coord = GRAPH_X_MID, self.GRAPH_Y-TITLE_SPACE, 0
        self._plotter.text(
            self._image, title, coord, 
            TITLE_DEPTH, TITLE_FONTSIZE, TITLE_COLOR,
            center_x=True, wrong_y=True, min_y=False)
            
    def draw_labels(
        self, xlabel, ylabel, zlabel, draw_x_axis, draw_y_axis, draw_z_axis,
        xtick_at, ytick_at, ztick_at, label_size=1.0):
        import graphconst as gc
        
        # Labels sit on the background.
        LABEL_DEPTH = 0.1*self.UNIT
        # Color (0.6, 0.6, 0.6) is too light for pybinreg predictions.
        LABEL_COLOR = 0.2, 0.2, 0.2
        LABEL_FONTSIZE = 2.5*self.UNIT*label_size
        LABEL_SPACE = 2.0*self.UNIT   # space between axis and label
        
        GRAPH_X_MID = self.GRAPH_X+self.GRAPH_WIDTH/2.0
        GRAPH_Y_MID = self.GRAPH_Y+self.GRAPH_HEIGHT/2.0
        GRAPH_Z_MID = self.GRAPH_Z+self.GRAPH_DEPTH/2.0

        if xlabel:
            max_y = True
            coord = [GRAPH_X_MID, LABEL_SPACE, 0]
            if not xtick_at&gc.BACK:
                coord[2] = self.GRAPH_Z + self.GRAPH_DEPTH - LABEL_DEPTH/2.0
            if not draw_x_axis:
                coord[1] = self.GRAPH_Y + self.GRAPH_HEIGHT + LABEL_SPACE
                max_y = False
            self._plotter.text(
                self._image, xlabel, coord, 
                LABEL_DEPTH, LABEL_FONTSIZE, LABEL_COLOR,
                center_x=True, max_y=max_y)
        if ylabel:
            min_x = True
            coord = [-LABEL_SPACE, GRAPH_Y_MID, 0]
            if not ytick_at&gc.BACK:
                coord[2] = self.GRAPH_Z + self.GRAPH_DEPTH - LABEL_DEPTH/2
            if not draw_y_axis:
                coord[0] = self.GRAPH_X - LABEL_SPACE
                min_x = False
            self._plotter.text90(
                self._image, ylabel, coord, 
                LABEL_DEPTH, LABEL_FONTSIZE, LABEL_COLOR,
                center_y=True, wrong_x=True, min_x=min_x)
        if zlabel and draw_z_axis:
            coord = [
                self.GRAPH_X, LABEL_SPACE, GRAPH_Z_MID-LABEL_DEPTH/2]
                #self.GRAPH_X, self.GRAPH_Y+self.GRAPH_HEIGHT+LABEL_SPACE,
                #GRAPH_Z_MID-LABEL_DEPTH/2]
            if not ztick_at&gc.LEFT:
                coord[0] = self.GRAPH_X + self.GRAPH_WIDTH + LABEL_SPACE
            num_to_skip = int(bool(xlabel)) + int(bool(ylabel))
            self._plotter.text(
                self._image, zlabel, coord, 
                LABEL_DEPTH, LABEL_FONTSIZE, LABEL_COLOR,
                center_x=True, max_y=num_to_skip)
                #center_x=True)

    def draw_points(
        self, points, color, shape, error_bars,
        onpoint_labels, onpoint_label_color, onpoint_label_size,
        overpoint_labels, overpoint_label_size, 
        shadow, point_size):
        # color should be parallel to points.
        # onpoint_label_color should be parallel to onpoint_labels.
        import colorlib
        import graphconst as gc

        assert len(points) == len(color)
        assert len(points) == len(shape)
        assert len(onpoint_labels) == len(onpoint_label_color)
        assert len(point_size) == len(points)
        
        #RADIUS = 0.8*self.UNIT*point_size  # Radius of each point.
        RADIUS = 0.8*self.UNIT             # Radius of each point.
        HEIGHT = 0                         # Bottom of each point.
        if self._zstack:
            HEIGHT = max([z+depth for (name, z, depth) in self._zstack])
        
        ERROR_RADIUS = 0.2*RADIUS
        ERROR_WIDTH = 1.1*2*RADIUS

        #ON_LABEL_FONTSIZE = RADIUS*1.25*onpoint_label_size
        #ON_LABEL_FONTSIZE = RADIUS*0.8*onpoint_label_size
        ON_LABEL_FONTSIZE = RADIUS*1.0*onpoint_label_size
        #ON_LABEL_DEPTH = RADIUS/2.0
        ON_LABEL_DEPTH = 0.1*self.UNIT*onpoint_label_size
        ON_LABEL_COLOR = 1, 1, 1

        OVER_LABEL_FONTSIZE = 1.35*RADIUS*overpoint_label_size
        OVER_LABEL_ALTFONTSIZE = OVER_LABEL_FONTSIZE * 0.75
        #OVER_LABEL_DEPTH = RADIUS/2.0
        OVER_LABEL_DEPTH = RADIUS*0.25*overpoint_label_size
        OVER_LABEL_COLOR = 0.2, 0.2, 0.2
        #OVER_LABEL_SPACE = 0.50 * RADIUS   # between edge of point and label
        OVER_LABEL_SPACE = OVER_LABEL_FONTSIZE*0.25  # btwn point edge & label

        point_coords = [None] * len(points)
        for i in range(len(points)):
            x_virt, y_virt = points[i][:2]
            z_virt = None
            if len(points[i]) == 3:
                z_virt = points[i][2]

            RAD = RADIUS*point_size[i]
            x_pix, y_pix = self._virtx2pix(x_virt), self._virty2pix(y_virt)
            z_pix = HEIGHT   # Bottom of each point.
            if z_virt is not None:
                z_pix = self._virtz2pix(z_virt)-RAD
            point_coords[i] = (
                x_pix-RAD, y_pix-RAD, z_pix, RAD*2, RAD*2, RAD*2)

            # Don't draw points that are off the drawing area.
            if x_pix<self.GRAPH_X or x_pix >= self.GRAPH_X+self.GRAPH_WIDTH:
                continue
            if y_pix<self.GRAPH_Y or y_pix >= self.GRAPH_Y+self.GRAPH_HEIGHT:
                continue

            col = color[i] or (0, 0, 0)
            on_label_col = onpoint_label_color[i]
            if on_label_col is None:
                #on_label_col = colorlib.choose_contrasting_color(col)
                on_label_col = colorlib.choose_contrasting_bw(col)

            # Draw the point.
            coord = x_pix, y_pix, z_pix+RAD
            self._plotter.sphere(
                self._image, coord, RAD, col, shape=shape[i],
                finish=gc.METALLIC, shadow=shadow)
            x = "points", z_pix, RAD*2
            if x not in self._zstack:
                self._zstack.append(x)

            # Draw the error bars.
            err_l = err_u = err_total = 0
            if error_bars:
                virt_height = self._virt_ylim[1] - self._virt_ylim[0]
                y_scale = self.GRAPH_HEIGHT / float(virt_height)
                x1, x2 = error_bars[i]
                err_l, err_u = x1*y_scale, x2*y_scale
                err_total = err_l + err_u
            if err_total:
                # Make sure the error bars don't go out of the plot.
                y_min = max(y_pix-err_u, self.GRAPH_Y)
                y_max = min(y_pix+err_l, self.GRAPH_Y+self.GRAPH_HEIGHT)
                err_height = y_max - y_min

                ER = ERROR_RADIUS * point_size[i]
                EW = ERROR_WIDTH * point_size[i]
                # Vertical error bar, from top to bottom.
                self._plotter.cylinder(
                    self._image, (x_pix, y_min, z_pix+RAD),
                    (0, err_height, 0), ER, col,
                    finish=gc.METALLIC)
                # Top horizontal bar.
                if y_pix-err_u >= self.GRAPH_Y:
                    self._plotter.cylinder(
                        self._image,
                        (x_pix-EW/2, y_pix-err_u, z_pix+RAD),
                        (EW, 0, 0), ER, col,
                        finish=gc.METALLIC)
                # Bottom horizontal bar.
                if y_pix+err_l < self.GRAPH_Y+self.GRAPH_HEIGHT:
                    self._plotter.cylinder(
                        self._image,
                        (x_pix-EW/2, y_pix+err_l, z_pix+RAD),
                        (EW, 0, 0), ER, col,
                        finish=gc.METALLIC)

                # Include the error bars into point_coords.
                radius = max(RAD, ER)
                y = min(y_pix-err_u-ER, y_pix-RAD)
                h = max(err_total+ER*2, RAD*2)
                point_coords[i] = (
                    x_pix-radius, y, z_pix, radius*2, h, radius*2)

            # Draw the labels on the points.
            label = onpoint_labels[i]
            if label:
                self._plotter.text(
                    self._image, label, 
                    (x_pix, y_pix, z_pix+RAD*2),
                    ON_LABEL_DEPTH, ON_LABEL_FONTSIZE, on_label_col,
                    center_x=True, center_y=True)

        # Position the labels around the points.
        plot_coords = (
            self.GRAPH_X, self.GRAPH_Y, self.GRAPH_Z,
            self.GRAPH_WIDTH, self.GRAPH_HEIGHT, self.GRAPH_DEPTH)
        label_sizes = [(0, 0, 0)] * len(points)
        altlabel_sizes = [(0, 0, 0)] * len(points)
        for i in range(len(points)):
            label = overpoint_labels[i]
            if not label:
                continue
            lw, lh = self._plotter.get_text_size(label, OVER_LABEL_FONTSIZE)
            ld = OVER_LABEL_DEPTH
            label_sizes[i] = lw, lh, ld
            lw, lh = self._plotter.get_text_size(label, OVER_LABEL_ALTFONTSIZE)
            altlabel_sizes[i] = lw, lh, ld
        margin = OVER_LABEL_SPACE
        x = position_labels(
            plot_coords, point_coords, label_sizes, altlabel_sizes, margin)
        label_positions, label_orientations, label_alternates, x = x

        # Draw the labels around the points.
        for i in range(len(points)):
            label = overpoint_labels[i]
            if not label:
                continue
            px, py, pz, pw, ph, pd = point_coords[i]
            #draw_above = y_pix > self.GRAPH_Y + self.GRAPH_HEIGHT/2

            lz = pz + pd + OVER_LABEL_DEPTH
            center_x = wrong_x = False
            center_y = wrong_y = False
            if label_positions[i] == gc.LABEL_TOP:
                lx = px + pw/2.0
                ly = py - OVER_LABEL_SPACE
                center_x = True
                wrong_y = True  # coordinate should be the bottom.
            elif label_positions[i] == gc.LABEL_BOTTOM:
                lx = px + pw/2.0
                ly = py + ph + OVER_LABEL_SPACE
                center_x = True
            elif label_positions[i] == gc.LABEL_LEFT:
                lx = px - OVER_LABEL_SPACE
                ly = py + ph/2.0
                wrong_x = True
                center_y = True
            elif label_positions[i] == gc.LABEL_RIGHT:
                lx = px + pw + OVER_LABEL_SPACE
                ly = py + ph/2.0
                center_y = True
            elif label_positions[i] == gc.LABEL_TOP_LEFT:
                lx = px - 0*OVER_LABEL_SPACE
                ly = py - 0*OVER_LABEL_SPACE
                wrong_x = True
                wrong_y = True
            elif label_positions[i] == gc.LABEL_TOP_RIGHT:
                lx = px + pw + 0*OVER_LABEL_SPACE
                ly = py - 0*OVER_LABEL_SPACE
                wrong_y = True
            elif label_positions[i] == gc.LABEL_BOTTOM_LEFT:
                lx = px - 0*OVER_LABEL_SPACE
                ly = py + ph + 0*OVER_LABEL_SPACE
                wrong_x = True
            elif label_positions[i] == gc.LABEL_BOTTOM_RIGHT:
                lx = px + pw + 0*OVER_LABEL_SPACE
                ly = py + ph + 0*OVER_LABEL_SPACE
            else:
                raise AssertionError, "Unknown label position"
            #coord = x_pix, label_y, z_pix+RADIUS*2+OVER_LABEL_DEPTH/2.0
            coord = lx, ly, lz
            fn = self._plotter.text
            #if overpoint_label_orientation[i] == gc.LABEL_VERTICAL:
            if label_orientations[i] == gc.LABEL_VERTICAL:
                fn = self._plotter.text90
            fontsize = OVER_LABEL_FONTSIZE
            if label_alternates[i]:
                fontsize = OVER_LABEL_ALTFONTSIZE
            fn(
                self._image, label, coord,
                OVER_LABEL_DEPTH, fontsize, OVER_LABEL_COLOR,
                center_x=center_x, wrong_x=wrong_x,
                center_y=center_y, wrong_y=wrong_y,
                background=True)
            
    def draw_line(self, points, color, shadow, same_height, line_size):
        import math
        import graphconst as gc

        RADIUS = 0.40 * self.UNIT * line_size
        HEIGHT = 0.0     # Bottom of each line.
        if self._zstack:
            if same_height:
                # Make same height as highest thing.
                HEIGHT = max([z for (name, z, depth) in self._zstack])
            else:
                # Make higher than highest thing.
                HEIGHT = max([z+depth for (name, z, depth) in self._zstack])
                
        if color is None:
            color = (0, 0, 0)
        
        # For convenience.
        GRAPH_X, GRAPH_Y, GRAPH_Z = self.GRAPH_X, self.GRAPH_Y, self.GRAPH_Z
        GRAPH_WIDTH, GRAPH_HEIGHT = self.GRAPH_WIDTH, self.GRAPH_HEIGHT
        GRAPH_DEPTH = self.GRAPH_DEPTH

        for i in range(len(points)-1):
            xv1, yv1 = points[i][:2]
            xv2, yv2 = points[i+1][:2]
            zv1 = zv2 = None
            if len(points[i]) == 3:
                zv1 = points[i][2]
            if len(points[i+1]) == 3:
                zv2 = points[i+1][2]

            xp1, yp1 = self._virtx2pix(xv1), self._virty2pix(yv1)
            xp2, yp2 = self._virtx2pix(xv2), self._virty2pix(yv2)
            zp1 = zp2 = HEIGHT+RADIUS   # Center of each line.
            if zv1 is not None:
                zp1 = self._virtz2pix(zv1)
            if zv2 is not None:
                zp2 = self._virtz2pix(zv2)

            # Make sure all the points are within the drawing area.
            # Both X-coordinates are off the plot.
            if max(xp1, xp2) < GRAPH_X or min(xp1, xp2) >= GRAPH_X+GRAPH_WIDTH:
                continue
            # Both Y-coordinates are off the plot.
            if max(yp1, yp2) < GRAPH_Y or min(yp1, yp2) >=GRAPH_Y+GRAPH_HEIGHT:
                continue
            # BUG: Check Z-coord.
            
            # At this point, at least one of the points are within the
            # drawing area.
            EPS = 1E-8
            if abs(xp2-xp1) < EPS and abs(yp2-yp1) < EPS:
                # Just a point.  Ignore.
                continue
            elif abs(xp2-xp1) < EPS:
                # Slope is infinity, straight line up and down.
                if xp1 < GRAPH_X or xp1 > GRAPH_X+GRAPH_WIDTH:
                    continue
                yp1 = min(max(yp1, GRAPH_Y), GRAPH_Y+GRAPH_HEIGHT)
                yp2 = min(max(yp2, GRAPH_Y), GRAPH_Y+GRAPH_HEIGHT)
            elif abs(yp2-yp1) < EPS:
                # Slope is 0, straight line left and right.
                if yp1 < GRAPH_Y or yp1 > GRAPH_Y+GRAPH_HEIGHT:
                    continue
                xp1 = min(max(xp1, GRAPH_X), GRAPH_X+GRAPH_WIDTH)
                xp2 = min(max(xp2, GRAPH_X), GRAPH_X+GRAPH_WIDTH)
            else:
                m = float((yp2-yp1)) / (xp2-xp1)
                b = yp1 - m * xp1
                # Make sure the x-coords are in the plot.
                min_x, max_x = GRAPH_X, GRAPH_X+GRAPH_WIDTH
                min_y, max_y = m*min_x + b, m*max_x + b
                if xp1 < min_x:
                    xp1, yp1 = min_x, min_y
                if xp1 > max_x:
                    xp1, yp1 = max_x, max_y
                if xp2 < min_x:
                    xp2, yp2 = min_x, min_y
                if xp2 > max_x:
                    xp2, yp2 = max_x, max_y
                # Make sure the y-coords are in the plot.
                min_y, max_y = GRAPH_Y, GRAPH_Y+GRAPH_HEIGHT
                min_x, max_x = (min_y-b)/m, (max_y-b)/m
                if yp1 < min_y:
                    xp1, yp1 = min_x, min_y
                if yp1 > max_y:
                    xp1, yp1 = max_x, max_y
                if yp2 < min_y:
                    xp2, yp2 = min_x, min_y
                if yp2 > max_y:
                    xp2, yp2 = max_x, max_y

            #x = xv1, yv1, xv2, yv2, xp1, yp1, xp2, yp2
            #print "\t".join(map(str, x))
            if abs(xp2-xp1) < EPS and abs(yp2-yp1) < EPS:
                # Just a point.  Ignore.
                continue

            # Draw a line.
            linelen = math.sqrt((xp2-xp1)**2 + (yp2-yp1)**2)
            angle = 90
            if abs(xp2-xp1) > 1E-8:
                angle = math.degrees(math.atan(float(yp2-yp1)/(xp2-xp1)))
            coord = xp1, yp1, zp1
            extent = xp2-xp1, yp2-yp1, zp2-zp1
            #if i == 0:
            #    self._plotter.sphere(
            #        self._image, coord, RADIUS, color,
            #        finish=gc.METALLIC, shadow=shadow)
            #print "PLOTTING", coord, extent
            self._plotter.cylinder(
                self._image, coord, extent, RADIUS, color, 
                finish=gc.METALLIC, shadow=shadow, blob=True)
            #coord = xp2, yp2, zp2
            #self._plotter.sphere(
            #    self._image, coord, RADIUS, color, 
            #    finish=gc.METALLIC, shadow=shadow)

            x = "line", max(zp1, zp2)-RADIUS, RADIUS*2
            if x not in self._zstack:
                self._zstack.append(x)

    def draw_bars(self, bars, barwidth, color, shadow,
                  overbar_label, overbar_label_size,
                  underbar_label, underbar_label_size,
                  width_size=1.0):
        import math
        import graphconst as gc

        assert len(overbar_label) == len(bars)

        DEPTH = 0.40 * self.UNIT
        HEIGHT = 0.0
        zstack = [x for x in self._zstack if x[0] != "bars"]
        if zstack:
            HEIGHT = max([z+depth for (name, z, depth) in zstack])

        LABEL_COLOR = 0.4, 0.4, 0.4
        OVERBAR_LABEL_FONTSIZE = 1.0*self.UNIT*overbar_label_size
        UNDERBAR_LABEL_FONTSIZE = 1.0*self.UNIT*overbar_label_size
        LABEL_DEPTH = 0.1*self.UNIT
        LABEL_SPACE = 0.1*self.UNIT   # space between bar and label

        if color is None:
            color = 0, 0, 0

        # For convenience.
        GRAPH_X, GRAPH_Y, GRAPH_Z = self.GRAPH_X, self.GRAPH_Y, self.GRAPH_Z
        GRAPH_WIDTH, GRAPH_HEIGHT = self.GRAPH_WIDTH, self.GRAPH_HEIGHT
        GRAPH_DEPTH = self.GRAPH_DEPTH

        for i in range(len(bars)):
            assert len(bars[i]) in [2, 3, 4]
            zc = None
            yminc = 0.0
            if len(bars[i]) == 2:
                xc, ymaxc = bars[i]
            elif len(bars[i]) == 3:
                #xc, yminc, ymaxc = bars[i]
                xc, ymaxc, zc = bars[i]
            else:
                #xc, yminc, ymaxc, zc = bars[i]
                xc, ymaxc, zc, yminc = bars[i]
            #assert yminc <= ymaxc, "%g %g" % (yminc, ymaxc)
            if yminc > ymaxc:
                yminc, ymaxc = ymaxc, yminc
            xc -= barwidth/2.0

            xp = self._virtx2pix(xc)
            yp = self._virty2pix(ymaxc)
            ymaxp = self._virty2pix(yminc)
            zp = HEIGHT
            if self._virt_zlim is not None:
                zp = self._virtz2pix(zc)

            xlen = self._virt_xlim[1] - self._virt_xlim[0]
            width = barwidth / float(xlen) * GRAPH_WIDTH * width_size
            height = ymaxp - yp

            # Make sure the points are within the drawing area.
            if xp < GRAPH_X:
                width = width - (GRAPH_X-xp)
                xp = GRAPH_X
            if xp+width > GRAPH_X+GRAPH_WIDTH:
                width = width - ((xp+width)-(GRAPH_X+GRAPH_WIDTH))
            if yp < GRAPH_Y:
                height = height - (GRAPH_Y-yp)
                yp = GRAPH_Y
            if yp+height > GRAPH_Y+GRAPH_HEIGHT:
                height = height - ((yp+height)-(GRAPH_Y+GRAPH_HEIGHT))
            if width <= 0 or height <= 0:
                continue
            
            # Draw the box.
            coord = xp, yp, zp
            extent = width, height, DEPTH
            #finish = gc.METALLIC
            finish = gc.ROUGH
            self._plotter.box(
                self._image, coord, extent, color, finish=finish,shadow=shadow)

            # Draw the label.
            label = overbar_label[i]
            if label:
                yc = yp - LABEL_SPACE
                wrong_y = True
                if yminc < 0:  # bar is negative
                    yc = yp + height + LABEL_SPACE
                    wrong_y = False
                coord = xp+width/2.0, yc, zp
                self._plotter.text(
                    self._image, label, coord,
                    LABEL_DEPTH, OVERBAR_LABEL_FONTSIZE, LABEL_COLOR,
                    center_x=True, wrong_y=wrong_y, vertical=True)
            label = underbar_label[i]
            if label:
                yc = yp + height + LABEL_SPACE
                wrong_y = False
                if yminc < 0:  # bar is negative
                    yc = yp - LABEL_SPACE
                    wrong_y = True
                coord = xp+width/2.0, yc, zp
                self._plotter.text(
                    self._image, label, coord,
                    LABEL_DEPTH, UNDERBAR_LABEL_FONTSIZE, LABEL_COLOR,
                    center_x=True, wrong_y=wrong_y, vertical=True)

            x = "bars", zp, DEPTH
            if x not in self._zstack:
                self._zstack.append(x)

    def write(self, handle, **keywds):
        if type(handle) is type(""):
            handle = open(handle, 'w')
        return self._plotter.write(self._image, handle, **keywds)
    
    def _virtx2pix(self, virt_x):
        # The user provides virtual coordinates.
        # Convert those coordinates to pixels on the plot.
        xlim = self._virt_xlim
        virt_min, virt_len = xlim[0], xlim[1]-xlim[0]
        pix_min, pix_len = self.GRAPH_X, self.GRAPH_WIDTH
        pix_x = float(virt_x-virt_min)/virt_len*pix_len + pix_min
        return pix_x

    def _virty2pix(self, virt_y):
        # Virtual coordinates have a Y that increases upwards.  Our
        # coordinates have a Y that increases downwards.
        ylim = self._virt_ylim
        virt_min, virt_len = ylim[0], ylim[1]-ylim[0]
        pix_min, pix_len = self.GRAPH_Y+self.GRAPH_HEIGHT, self.GRAPH_HEIGHT
        pix_y = pix_min - float(virt_y-virt_min)/virt_len*pix_len
        return pix_y

    def _virtz2pix(self, virt_z):
        # The user provides virtual coordinates.
        # Convert those coordinates to pixels on the plot.
        zlim = self._virt_zlim
        virt_min, virt_len = zlim[0], zlim[1]-zlim[0]
        pix_min, pix_len = self.GRAPH_Z, self.GRAPH_DEPTH
        pix_z = float(virt_z-virt_min)/virt_len*pix_len + pix_min
        return pix_z

def scatter(*args, **keywds):
    """Return the Graph object.
    
    points           List of (x, y) coordinates for the points.
                     Also can be (x, y, z).
    color            List of the color for each point (<r>,<g>,<b>) from 0-1.
    shape            List of DEFAULT, CIRCLE, SQUARE, or DIAMOND.
    shadow
    onpoint_label    What label to put on each point.
    onpoint_label_color
    onpoint_label_size
    overpoint_label  What label to put over each point.
    overpoint_label_size
    error_bar        List of error or (errlow, errhi).
                     The errors are the heights of the error bars.
    point_size       Scales the size of each point.
    More from _make_graph.
    
    """
    #from StringIO import StringIO
    #import graphconst as gc

    assert len(args) == 1, "Specify points"
    points, = args
    color = keywds.get("color", None)
    shape = keywds.get("shape", None)
    shadow = keywds.get("shadow", None)
    onpoint_label = keywds.get("onpoint_label", None)
    onpoint_label_color = keywds.get("onpoint_label_color", None)
    onpoint_label_size = keywds.get("onpoint_label_size", 1.0)
    overpoint_label = keywds.get("overpoint_label", None)
    overpoint_label_size = keywds.get("overpoint_label_size", 1.0)
    #overpoint_label_orientation = keywds.get(
    #    "overpoint_label_orientation", gc.LABEL_HORIZONTAL)
    error_bar = keywds.get("error_bar", None)
    point_size = keywds.get("point_size", 1.0)
    graph = keywds.get("graph", None)
    plotter = keywds.get("plotter", None)

    assert points, "No points provided for scatter plot."

    # Check the inputs
    point_size = _set_default_point_size(point_size, len(points))
    color = _set_default_color(color, len(points))
    shape = _set_default_shape(shape, len(points))
    onpoint_label = _set_default_label(onpoint_label, len(points))
    onpoint_label_color = _set_default_color(
        onpoint_label_color, len(onpoint_label))
    overpoint_label = _set_default_label(overpoint_label, len(points))
    #overpoint_label_orientation = _set_default_label(
    #    overpoint_label_orientation, len(points))
    error_bar = _set_default_error_bar(error_bar)
    assert len(point_size) == len(points)
    assert len(color) == len(points)
    assert len(shape) == len(points)
    assert len(onpoint_label) == len(points)
    assert len(onpoint_label_color) == len(onpoint_label)
    assert len(overpoint_label) == len(points)
    assert not error_bar or len(error_bar) == len(points)

    X = [x[0] for x in points]
    Y = [x[1] for x in points]
    Z = []
    for x in points:
        if len(x) == 3:
            Z.append(x[2])
    if graph is None:
        Z_min = Z_max = None
        if Z:
            Z_min, Z_max = min(Z), max(Z)
        graph = _make_graph(
            plotter, min(X), max(X), min(Y), max(Y), Z_min, Z_max, **keywds)
    graph.draw_points(
        points, color, shape, error_bar,
        onpoint_label, onpoint_label_color, onpoint_label_size, 
        overpoint_label, overpoint_label_size, 
        shadow, point_size)
    return graph

def line(*args, **keywds):
    """Return the Graph object.
    
    line1        List of the (x, y) coordinates of the points.
    line2        ...
    color        List of colors, 1 for each line (<r>,<g>,<b>) from 0-1.
    shape        List of DEFAULT, CIRCLE, SQUARE, or DIAMOND.
    shadow
    line_size    Scales the size of the line.
    draw_points  Whether to draw the points or not.
    point_size   Scales the size of each point.
    same_height  Draw all lines at the same height (default False).
    
    """
    import itertools

    #KNOWN_KEYWDS = [
    #    "color", "shape", "draw_points", "point_size",
    #    "xlim", "ylim",
    #    "xtick", "xtick_label", "ytick", "ytick_label", "tick_size",
    #    "xlabel", "ylabel", "label_size", "font",
    #    "width", "height"
    #    ]
    #for word in keywds:
    #    assert word in KNOWN_KEYWDS, "I don't know the keyword: %s" % word
    assert args, "No lines given."

    color = keywds.get("color", None)
    shape = keywds.get("shape", None)
    shadow = keywds.get("shadow", None)
    line_size = keywds.get("line_size", 1.0)
    draw_points = keywds.get("draw_points", None)
    point_size = keywds.get("point_size", 1.0)
    same_height = keywds.get("same_height", False)
    graph = keywds.get("graph", None)
    plotter = keywds.get("plotter", None)

    # Check the inputs
    lines = args
    for line in lines:
        for x in line:
            assert len(x) in [2, 3], str(x)
    X = list(itertools.chain.from_iterable(
        [[x[0] for x in line] for line in lines]))
    Y = list(itertools.chain.from_iterable(
        [[x[1] for x in line] for line in lines]))
    Z = []
    for line in lines:
        for x in line:
            if len(x) == 3:
                Z.append(x[2])
    
    point_size = _set_default_point_size(point_size, len(lines))
    color = _set_default_color(color, len(lines))
    shape = _set_default_shape(shape, len(lines))
    assert len(point_size) == len(lines)
    assert len(color) == len(lines)
    assert len(shape) == len(lines)

    # The points on the line should be smaller than the points in the
    # scatter plot.
    point_size = [x*0.80*line_size for x in point_size]

    if not graph:
        Z_min = Z_max = None
        if Z:
            Z_min, Z_max = min(Z), max(Z)
        graph = _make_graph(
            plotter, min(X), max(X), min(Y), max(Y), Z_min, Z_max, **keywds)

    for i, line in enumerate(lines):
        graph.draw_line(line, color[i], shadow, same_height, line_size)
        if not draw_points:
            continue
        # Should let the user specify which points to plot.
        col = [color[i]] * len(line)
        shp = [shape[i]] * len(line)
        default = [None] * len(line)
        graph.draw_points(
            line, col, shp,
            None,                    # error_bars
            default, default, 1.0,   # onpoint_labels, color, size
            default, 1.0,            # overpoint_labels, size
            shadow, point_size)
    return graph

def bar(*args, **keywds):
    """Return the Graph object.

    series1        List of the (x, y) coordinates of this set of bars.
                   Can also be (x, y, z).
                   Can also be (x, y, z, y_min).
    series2        ...
    color          List of colors, 1 for each bar (<r>,<g>,<b>) from 0-1.
    shadow         Whether to draw a shadow.
    width_size     Width of bar.  1.0 means bars do not overlap.
    overbar_label  Label to draw over each bar.
    overbar_label_size
    underbar_label  Label to draw under each bar.
    underbar_label_size
    bar_label      Label for each bar.
    bar_label_size

    """
    # TODO:
    # ERROR BAR
    # Label multiple series.
    # detect whether xtick_label should be vertical
    color = keywds.get("color", None)
    shadow = keywds.get("shadow", None)
    width_size = keywds.get("width_size", 1.0)
    graph = keywds.get("graph", None)
    plotter = keywds.get("plotter", None)
    bar_label = keywds.get("bar_label", None)
    bar_label_size = keywds.get("bar_label_size", 1.0)
    overbar_label = keywds.get("overbar_label", None)
    overbar_label_size = keywds.get("overbar_label_size", 1.0)
    underbar_label = keywds.get("underbar_label", None)
    underbar_label_size = keywds.get("underbar_label_size", 1.0)

    # Check the inputs
    assert args, "No bars given."
    all_series = args
    X, Y, Z = [], [], []
    for series in all_series:
        for bar in series:
            assert len(bar) in [2, 3, 4]
            X.append(bar[0])
            Y.append(bar[1])
            if len(bar) == 2:
                Y.append(0.0)
            elif len(bar) == 3:
                Y.append(bar[2])
            elif len(bar) == 4:
                Y.append(bar[2])
                Z.append(bar[3])
                
    color = _set_default_color(color, len(all_series))
    assert len(color) == len(all_series)

    # BUG: label only applies to the first series.
    assert "xtick" not in keywds
    assert "xtick_label" not in keywds
    assert "vertical_xtick_label" not in keywds
    xtick = xtick_label = None
    if bar_label:
        series = all_series[0]
        xtick = [bar[0] for bar in series]
        xtick_label = _set_default_label(bar_label, len(all_series[0]))
    overbar_label = _set_default_label(overbar_label, len(all_series[0]))
    underbar_label = _set_default_label(underbar_label, len(all_series[0]))

    barwidth = _set_bar_width(X)

    if graph is None:
        X_min = min(X) - barwidth/2.0
        X_max = max(X) + barwidth/2.0
        Z_min = Z_max = None
        if Z:
            Z_min, Z_max = min(Z), max(Z)
        # BUG: label_size parameter is ignored.
        graph = _make_graph(
            plotter, X_min, X_max, min(Y), max(Y), Z_min, Z_max,
            xtick=xtick, xtick_label=xtick_label, label_size=bar_label_size,
            vertical_xtick_label=True, **keywds)
        
    for i, series in enumerate(all_series):
        graph.draw_bars(
            series, barwidth, color[i], shadow,
            overbar_label, overbar_label_size,
            underbar_label, underbar_label_size,
            width_size=width_size)

    return graph

def histogram(X, *args, **keywds):
    """Return a Graph object.

    X            List of the values to plot.
    breaks       Number of breaks
    frequency    True/False.  Plot the frequency or the density.
    bar_size     1.0 means bars do not overlap.

    """
    # TODO: Should recognize if X is all integers and not allow any
    # bars at sub-integer frequencies.
    breaks = keywds.get("breaks", None)
    frequency = keywds.get("frequency", True)
    bar_size = keywds.get("bar_size", 1.0)
    
    x = calc_histogram(X, breaks=breaks)
    breaks, counts, density = x

    series = counts
    if not frequency:
        series = density
    # Convert to (x, y) coordinates.
    series = zip(breaks[:-1], series)
        
    graph = bar(series, width_size=bar_size, **keywds)
    return graph

def place_breaks(v_min, v_max, num_breaks=10, delta=None):
    assert v_min < v_max, "%g %g" % (v_min, v_max)
    assert num_breaks > 0 and num_breaks <= 1000

    breaks = place_ticks(v_min, v_max, num_ticks=num_breaks, delta=delta)
    assert len(breaks) >= 2
    delta = breaks[1] - breaks[0]
    if breaks[0] > v_min:
        breaks = [breaks[0]-delta] + breaks
    if breaks[-1] < v_max:
        breaks = breaks + [breaks[-1]+delta]
    assert breaks[0] >= v_min
    assert breaks[-1] <= v_max
    return breaks

def calc_histogram(X, breaks=None):
    import bisect
    
    if breaks is None:
        breaks = place_breaks(min(X), max(X))
    elif type(breaks) is type(0):
        breaks = place_breaks(min(X), max(X), num_breaks=breaks)

    counts = [0]*(len(breaks)-1)
    for x in X:
        # Find the correct bin.
        i = bisect.bisect_left(breaks, x)
        # Drop the items that are too big.  If x is exactly
        # breaks[-1], then will not count.
        if i >= len(counts):  # if x is too big.
            continue
        # Drop the items that are too small.
        if breaks[i] > x:     # if x is too small.
            continue
        # Increase the right count.
        counts[i] += 1
    sum_counts = float(sum(counts))
    density = [x/sum_counts for x in counts]
    return breaks, counts, density

def plot_heatmap(infile, outfile, xpix, ypix, **keywds):
    import os
    import sys
    import subprocess

    # Set the path to the same path as the outfile.
    path, x = os.path.split(outfile)
    cwd = None
    try:
        if path:
            cwd = os.getcwd()
            os.chdir(path)
    
        cmd = plot_heatmap_cmd(infile, outfile, xpix, ypix, **keywds)
        p = subprocess.Popen(
            cmd, shell=True, bufsize=0, stdin=subprocess.PIPE,
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, close_fds=True)
        w, r = p.stdin, p.stdout
        w.close()
        output = r.read()
    finally:
        if cwd:
            os.chdir(cwd)

    # Make sure the signature was generated correctly.  An error could
    # mean that arrayplot.py or cluster is missing.
    if not os.path.exists(outfile):
        print >>sys.stderr, output
        raise AssertionError, "Failed to make dataset.\n%s" % cmd
    return output

def plot_heatmap_cmd(
    infile, outfile, xpix, ypix, color=None,
    show_colorbar=None, show_grid=None, 
    scale=None, gain=None, no_autoscale=False, 
    gene_label=False, cluster_genes=False,
    gene_center=None, gene_normalize=None,
    array_label=False, cluster_arrays=False,
    python=None, arrayplot=None, cluster=None, libpath=None):
    import os
    import config

    python = python or config.python
    arrayplot = arrayplot or config.arrayplot or "arrayplot.py"
    # If arrayplot.py is not a real file, then it may need to be found
    # on the path.  In this case, make sure python is None.
    if not os.path.exists(arrayplot) and \
       not os.path.exists(os.path.realpath(arrayplot)):
        python = None
    color = color or "bild"

    cmd = [
        python,
        arrayplot,
        "-x %d" % xpix,
        "-y %d" % ypix,
        "-o %s" % outfile,
        "--color=%s" % color,
        ]
    if show_grid:
        cmd.append("--grid")
    if show_colorbar:
        cmd.append("--colorbar")
    if gene_label:
        cmd.append("--gl")
    if cluster_genes:
        cmd.append("-g")
    if gene_center:
        assert gene_center in ["mean", "median"]
        cmd.append("--gc=%s" % gene_center)
    if gene_normalize:
        assert gene_normalize in ["ss", "var"]
        cmd.append("--gn=%s" % gene_normalize)
    if array_label:
        cmd.append("--al")
    if cluster_arrays:
        cmd.append("-a")
    if scale:
        cmd.append("--scale=%s" % scale)
    if gain:
        cmd.append("--gain=%s" % gain)
    if no_autoscale:
        cmd.append("--no_autoscale")
    if cluster is not None:
        cmd.append("--cluster_app=%s" % cluster)
    if libpath:
        for path in libpath:
            cmd.append("--libpath=%s" % path)
    cmd.append("'%s'" % infile)

    # If python is None, then ignore it.
    cmd = [x for x in cmd if x is not None]
    cmd = " ".join(cmd)
    #print cmd
    return cmd

def find_tall_heatmap_size(
    nrow, ncol, min_box_height=None, min_box_width=None,
    max_box_height=None, max_box_width=None,
    max_total_height=None, max_total_width=None, height_width_ratio=None,
    max_megapixels=None):
    # Return tuple of the pixels for each box as (xpix, ypix).
    # Minimum sizes take precedence over maximum sizes.
    # max_megapixels take precedence over everything.
    import math

    # Set some defaults.
    min_box_height = min_box_height or 1.0
    min_box_width = min_box_width or 1.0
    max_total_height = max_total_height or 600
    max_total_width = max_total_width or 800
    max_box_height = max_box_height or max(min_box_height,
                                           float(max_total_height)/nrow)
    max_box_width = max_box_width or max(min_box_width,
                                         float(max_total_width)/ncol)
    height_width_ratio = height_width_ratio or 1.618  # HEIGHT/WIDTH

    min_box_height = float(min_box_height)
    min_box_width = float(min_box_width)
    max_box_height = float(max_box_height)
    max_box_width = float(max_box_width)
    max_total_height = float(max_total_height)
    max_total_width = float(max_total_width)
    height_width_ratio = float(height_width_ratio)
    
    # Start with the minimum size.
    ypix = min_box_height
    total_y = ypix * nrow
    # Use the height to find the right width.
    total_x = total_y / height_width_ratio
    xpix = float(total_x) / ncol
    #print total_x, total_y
    
    # If the width is too small, then set the width to the minimum and
    # use it to set the height.
    if xpix < min_box_width:
        xpix = min_box_width
        total_x = xpix * ncol
        total_y = total_x * height_width_ratio
        ypix = float(total_y) / nrow
        assert ypix >= min_box_height

    # Increase the size up to the maximum allowed.
    max_xpix = max(min(max_box_width, float(max_total_width) / ncol), 1.0)
    max_ypix = max(min(max_box_height, float(max_total_height) / nrow), 1.0)
    x_ratio = max_xpix / xpix
    y_ratio = max_ypix / ypix
    ratio = min(y_ratio, x_ratio)
    if ratio > 1.0:
        xpix, ypix = xpix*ratio, ypix*ratio
    #print xpix, ypix

    # Decrease the size to fit max_megapixels.
    megapixels = nrow * xpix * ncol * ypix / 1024.0 / 1024.0
    if max_megapixels is not None and megapixels > max_megapixels:
        assert max_megapixels > 0
        ratio = math.sqrt(max_megapixels / megapixels)
        invratio = 1.0 / ratio
        if invratio < xpix and invratio < ypix:
            # Try to just scale both dimensions.
            xpix = xpix * ratio
            ypix = ypix * ratio
        elif ncol*ypix > nrow*xpix:
            xpix = 1.0
            ypix = float(max_megapixels) / nrow / xpix / ncol * 1024.0 * 1024.0
        else:
            ypix = 1.0
            xpix = float(max_megapixels) / nrow / ncol / ypix * 1024.0 * 1024.0
        assert xpix >= 1 and ypix >= 1, "Plot exceeds megapixel limit."

    xpix, ypix = int(xpix), int(ypix)
    #height = nrow * ypix
    #width = ncol * xpix
    #print "DIM %dx%d (%d*%d, %d*%d)" % (height, width, nrow, ypix, ncol, xpix)
    return xpix, ypix

def find_wide_heatmap_size(
    nrow, ncol, min_box_height=None, min_box_width=None,
    max_box_height=None, max_box_width=None,
    max_total_height=None, max_total_width=None, height_width_ratio=None,
    max_megapixels=None):
    inv_height_width_ratio = height_width_ratio
    if inv_height_width_ratio is not None:
        inv_height_width_ratio = 1.0 / inv_height_width_ratio
    x = find_tall_heatmap_size(
        ncol, nrow,
        min_box_height=min_box_width, min_box_width=min_box_height,
        max_box_height=max_box_width, max_box_width=max_box_height,
        max_total_height=max_total_width, max_total_width=max_total_height,
        height_width_ratio=inv_height_width_ratio,
        max_megapixels=max_megapixels)
    xpix, ypix = x
    xpix, ypix = ypix, xpix
    #height = nrow * ypix
    #width = ncol * xpix
    #print "DIM %dx%d (%d*%d, %d*%d)" % (height, width, nrow, ypix, ncol, xpix)
    return xpix, ypix

def place_ticks(v_min, v_max, num_ticks=10, delta=None):
    import math
    
    if v_min == v_max:
        v_min = v_min * 0.95
        v_max = v_max * 1.05
    assert v_min < v_max, "%g %g" % (v_min, v_max)
    assert num_ticks > 0 and num_ticks <= 1000

    if delta is None:
        delta = _choose_tick_delta(v_min, v_max, num_ticks=num_ticks)

    # Do calculation in integers to 4 decimal places to avoid floating
    # point representation issues.
    x = math.log(delta, 10)
    multiplier = 10**(-int(x)+4)

    delta = int(delta * multiplier)
    tick_min, tick_max = int(v_min*multiplier), int(v_max*multiplier)

    # Round tick_min down to the nearest delta and tick_max up to the
    # nearest delta to cover the whole range.
    tick_min = tick_min - tick_min % delta
    if tick_max % delta != 0:
        tick_max = tick_max + (delta - tick_max % delta)

    ticks = [float(i)/multiplier for i in range(tick_min, tick_max+1, delta)]

    # Convert ticks to integers if possible.
    for t in ticks:
        if abs(float(t)-int(t)) >= 1E-50:
            break
    else:
        ticks = [int(x) for x in ticks]

    # Only keep the ticks that are within the range.
    ticks = [x for x in ticks if x >= v_min and x <= v_max]

    return ticks

def _make_graph(
    plotter, X_min, X_max, Y_min, Y_max, Z_min, Z_max,
    xlim=None, ylim=None, zlim=None, 
    xtick=None, xtick_label=True, ytick=None, ytick_label=True,
    ztick=None, ztick_label=True,
    vertical_xtick_label=False,
    xtick_at=None, ytick_at=None, ztick_at=None, tick_size=1.0,
    draw_x_axis=None, draw_y_axis=None, draw_z_axis=None,
    draw_box=None, 
    grid=None, title=None, title_size=1.0,
    xlabel=None, ylabel=None, zlabel=None, label_size=1.0,
    font=None, width=None, height=None, **keywds):
    # xlim         Tuple of (minimum x, maximum x).
    # xtick        Coordinates of the tickmarks on the x-axis (or True).
    # xtick_label  Labels for the tickmarks.  Should be parallel to xtick.
    # tick_size    Scales the size of the tick mark labels.
    #
    # grid         Whether to draw the grid in the background.
    # xlabel       Label for the X-axis.
    # ylabel
    # label_size   Scales the size of the labels.
    # font
    # 
    # width        Number of pixels wide for the plot.
    # height

    # keywds will absorb additional arguments...
    import graphconst as gc

    width = width or 1024
    height = height or 768
    if plotter is None:
        import povrayplot as plotter

    # Ignored!
    #if font:
    #    raise NotImplementedError

    # Check the inputs
    x = _set_default_axes(
        X_min, X_max, Y_min, Y_max, Z_min, Z_max, 
        xlim, ylim, zlim,
        xtick, ytick, ztick, xtick_label, ytick_label, ztick_label)
    (xlim, ylim, zlim,
     xtick, ytick, ztick, xtick_label, ytick_label, ztick_label) = x
     #print xlim, ylim, zlim

    if draw_x_axis is None:
        draw_x_axis = True
    if draw_y_axis is None:
        draw_y_axis = True
    if draw_z_axis is None and zlim:
        draw_z_axis = True

    if not draw_x_axis:
        xtick, xtick_label = [], None
    if not draw_y_axis:
        ytick, ytick_label = [], None
    if not draw_z_axis:
        ztick, ztick_label = [], None

    if ztick and grid is None:
        # None means no value set.  Set to False if really don't want grid.
        grid = True

    if xtick_at is None:
        xtick_at = gc.BACK         # default is bottom, back
        if ztick and grid:
            xtick_at = xtick_at^gc.BACK
    if ytick_at is None:
        ytick_at = gc.LEFT|gc.BACK # default is bottom, back
        if ztick and grid:
            ytick_at = ytick_at^gc.BACK
    if ztick_at is None:
        ztick_at = gc.LEFT         # default is bottom, left
        if grid:
            ztick_at = ztick_at^gc.LEFT

    #print draw_x_axis, draw_y_axis, draw_z_axis, draw_box
    graph = Graph(plotter, width, height, xlim, ylim, zlim)
    graph.draw_axes(draw_x_axis, draw_y_axis, draw_z_axis, draw_box)
    graph.draw_tick_marks(
        xtick, ytick, ztick, vertical_xtick_label,
        xtick_label, ytick_label, ztick_label,
        xtick_at, ytick_at, ztick_at, grid, tick_size=tick_size,
        label_size=label_size)
    graph.draw_title(title, title_size=title_size)
    graph.draw_labels(
        xlabel, ylabel, zlabel, 
        draw_x_axis, draw_y_axis, draw_z_axis,
        xtick_at, ytick_at, ztick_at, label_size=label_size)
    return graph

def _set_default_color(color, n):
    # Return an n-length list of colors.
    # 
    # As input, the possible values for color are:
    # (r, g, b)    Use this color for everything.  color from 0-1.
    # True         Use a default color scheme.
    # <function>   Use the function to generate a color scheme.
    # <false>      Everything is black (no color).
    import operator
    import colorlib

    if n == 0:
        return []
    assert n > 0
    if color is True:
        if n == 1:
            color = None
        else:
            color = colorlib.bild_colors(n)
    elif type(color) is type(_set_default_color):
        color = color(n)
    if not color:
        color = None
        
    # At this point, color can be:
    # 1.  None                                  None
    # 2.  (r, g, b)                             Color
    # 3.  [(r, g, b), ...]                      Sequence of colors.
    # 4.  [None, None, (r, g, b), None, ...]    Mix of None and colors.
    # 5.  [None, None, None, None]

    isnum = operator.isNumberType

    # Case 1: Convert None to list of black.
    if color is None:
        color = [(0, 0, 0)] * n
    # Case 2: Convert to list.
    if len(color) == 3 and \
           isnum(color[0]) and isnum(color[1]) and isnum(color[2]):
        color = [color] * n
    # Case 4, 5: Convert None to black.
    if operator.isSequenceType(color):
        for i in range(len(color)):
            if color[i] is None:
                color[i] = (0, 0, 0)

    # At this point, color should be list of (r, g, b).
    assert len(color) == n, "You specified %d colors, but I expected %d." % (
        len(color), n)
    for c in color:
        assert len(c) == 3, "Does not look like color: %s" % str(c)
        assert isnum(c[0]) and isnum(c[1]) and isnum(c[2]), \
               "Does not look like color: %s" % str(c)        
    return color

def _set_default_point_size(size, n):
    # Return a list of sizes.
    import operator
    if size is None:
        size = 1.0
    if type(size) is type(0) or type(size) is type(0.0):
        size = [size]
    assert operator.isSequenceType(size)
    if len(size) == 1:
        size = size * n
    assert len(size) == n
    return size

def _set_default_shape(shape, n):
    # Return a list of shapes.
    import operator
    import graphconst as gc
    if shape is None:
        shape = gc.DEFAULT
    if type(shape) is type(gc.DEFAULT):
        shape = [shape]
    assert operator.isSequenceType(shape)
    if len(shape) == 1:
        shape = shape * n
    assert len(shape) == n
    return shape

def _set_default_label(label, n):
    # Return a n-length list of labels.  As input, label can be:
    # List of values     None means don't label this point.
    # one value
    # None               Don't label anything.
    # []                 Don't label anything.
    import operator
    if label is None or label == []:
        label = [None]
    if not operator.isSequenceType(label):
        label = [label]
    assert operator.isSequenceType(label)
    if len(label) == 1:
        label = label * n
    assert len(label) == n
    return label

def _set_default_error_bar(error_bar):
    import operator
    
    if not error_bar:
        return None
    # Each error can be a single number or a tuple of (lower, upper).
    # Normalize so each of them are tuples of (lower, upper).
    error_bar = error_bar[:]
    for i in range(len(error_bar)):
        x = error_bar[i]
        if operator.isNumberType(x):
            err_l = err_u = x
        elif operator.isSequenceType(x):
            assert len(x) == 2
            err_l, err_u = x
        else:
            raise AssertionError, "Unknown error: %s" % repr(x)
        assert err_l >= 0 and err_u >= 0, "Errors: %s, %s" % (err_l, err_u)
        error_bar[i] = err_l, err_u
    return error_bar
    
def _set_default_axes(
    min_x, max_x, min_y, max_y, min_z, max_z, 
    xlim, ylim, zlim,
    xtick, ytick, ztick, xtick_label, ytick_label, ztick_label):
    # Initialize some variables.
    if xtick_label and not xtick:
        xtick = True
    if ytick_label and not ytick:
        ytick = True
    if ztick_label and not ztick:
        ztick = True

    # Set the tick marks.
    if xlim:
        assert len(xlim) == 2
        xtick = _set_default_tick(xtick, xlim[0], xlim[1])
    else:
        xtick = _set_default_tick(xtick, min_x, max_x)
    if ylim:
        assert len(ylim) == 2
        ytick = _set_default_tick(ytick, ylim[0], ylim[1])
    else:
        ytick = _set_default_tick(ytick, min_y, max_y)
    if zlim:
        assert len(zlim) == 2
        ztick = _set_default_tick(ztick, zlim[0], zlim[1])
    elif min_z is not None and max_z is not None:
        ztick = _set_default_tick(ztick, min_z, max_z)
    else:
        # No Z coordinates.
        ztick = None
        ztick_label = None

    # Set the tick mark labels.
    xtick_label = _set_default_tick_label(xtick_label, xtick)
    ytick_label = _set_default_tick_label(ytick_label, ytick)
    if ztick is not None:
        ztick_label = _set_default_tick_label(ztick_label, ztick)

    # Set the limits.
    if xtick:
        min_x, max_x = min(min_x, min(xtick)), max(max_x, max(xtick))
    if ytick:
        min_y, max_y = min(min_y, min(ytick)), max(max_y, max(ytick))
    if ztick:
        min_z, max_z = min(min_z, min(ztick)), max(max_z, max(ztick))
    xlim = _set_default_lim(xlim, min_x, max_x)
    ylim = _set_default_lim(ylim, min_y, max_y)
    if min_z is not None and max_z is not None:
        zlim = _set_default_lim(zlim, min_z, max_z)

    assert not xtick_label or len(xtick) == len(xtick_label)
    assert not ytick_label or len(ytick) == len(ytick_label)
    assert not ztick_label or len(ztick) == len(ztick_label)
    
    x = (xlim, ylim, zlim,
         xtick, ytick, ztick, xtick_label, ytick_label, ztick_label)
    return x

def _set_default_lim(lim, coord_min, coord_max):
    # Return (min, max).
    if lim:
        assert len(lim) == 2
        return lim
    coord_min, coord_max = float(coord_min), float(coord_max)
    x_len = coord_max - coord_min
    if x_len < 1E-99:  # max and min coord are the same.
        x_len = 1
    x_min = coord_min - x_len*0.05
    x_max = coord_max + x_len*0.05
    return x_min, x_max

def _set_default_tick(tick, coord_min, coord_max):
    import operator
    
    if not tick:
        tick = []
    elif operator.isSequenceType(tick):
        pass
    elif tick is True:
        tick = place_ticks(coord_min, coord_max, num_ticks=8)
    elif type(tick) is type(0):
        tick = place_ticks(coord_min, coord_max, num_ticks=tick)
    else:
        raise AssertionError, "I don't understand: %s" % str(tick)
    return tick

def _set_default_tick_label(label, tick):
    import operator
    if not tick:
        label = []
    elif operator.isSequenceType(label):
        assert len(label) == len(tick)
    elif label:
        # Figure out the number of digits after the decimal place.
        x = [max(len(str(abs(x)%1))-2, 0) for x in tick]
        digits = max(x)
        label = ["%.*f" % (digits, x) for x in tick]
    else:
        label = []
    return label

def _set_bar_width(X):
    # Width should be the minimum distance from one bar to the next.
    # Can draw the bars with no overlap.
    # This doesn't work for 3D plots since it doesn't take into
    # account the Z dimension.
    X = sorted(X)

    # Hack: Ignore bars that completely overlap.
    i = 0
    while i < len(X)-1:
        if abs(X[i]-X[i+1]) < 1E-50:
            del X[i]
        else:
            i += 1

    dists = [X[i]-X[i-1] for i in range(1, len(X))]
    width = min(dists)
    return width

def _choose_tick_delta(v_min, v_max, num_ticks=10):
    # Figure out the right delta to make at most num_ticks tick marks
    # between v_min and v_max.
    import math

    assert v_min < v_max
    assert num_ticks > 0 and num_ticks <= 100

    # Generate a list of the allowable DELTAs.  This will depend on
    # the range of the user's data.
    DELTAS = [0.5, 0.25, 0.20, 0.1]

    # Calculate the ideal delta, and choose the smallest DELTA that is
    # >= the ideal one.  Assume DELTAS sorted from largest to
    # smallest.  This means num_ticks is the maximum bound.
    range_ = v_max - v_min
    delta_ideal = float(range_) / num_ticks
    # delta_ideal   delta
    #   0.8           0.5
    #    12           10
    #    23           20
    #   0.02         0.02

    # Scale the DELTAs to be in the same range as delta_ideal.
    # delta_ideal should be smaller than max(DELTAS).    
    num_logs = int(math.ceil(math.log(delta_ideal / max(DELTAS), 10)))
    # Need to make DELTA this number of logs bigger.
    DELTAS = [x*10**num_logs for x in DELTAS]
    assert delta_ideal <= max(DELTAS)
    x = [x for x in DELTAS if x >= delta_ideal]
    delta = x[-1]

    return delta

def _overlaps_1d(x1, l1, x2, l2):
    if x1 > x2:
        x1, l1, x2, l2 = x2, l2, x1, l1
    return min(max(x1+l1-x2, 0), l2)

def _overlaps_3d(obj1, obj2):
    x1, y1, z1, w1, h1, d1 = obj1
    x2, y2, z2, w2, h2, d2 = obj2

    ox = _overlaps_1d(x1, w1, x2, w2)
    if not ox:  # optimization: short circuit the comparison.
        return 0
    oy = _overlaps_1d(y1, h1, y2, h2)
    if not oy:
        return 0
    oz = _overlaps_1d(z1, d1, z2, d2)
    return ox*oy*oz

def _overlaps(obj1, obj2):
    return _overlaps_3d(obj1, obj2)

def _perc_overlap(obj1, obj2):
    # Can rewrite this in C for a 10% speedup for pybinreg.py.
    x1, y1, z1, w1, h1, d1 = obj1
    x2, y2, z2, w2, h2, d2 = obj2
    area1 = w1*h1*d1
    area2 = w2*h2*d2
    x = _overlaps_3d(obj1, obj2)
    p = float(x) / min(area1, area2)
    EPS = 1E-10
    assert p >= 0.0-EPS and p <= 1.0+EPS, p
    return p

def _calc_label_coord(
    object, label, altlabel, position, orientation, alternate, margin):
    # object       (x, y, z, width, height, depth)
    # label        (width, height, depth)
    # altlabel     (width, height, depth)
    # position     from graphconst
    # orientation  from graphconst
    # alternate    Whether to use altlabel instead of label.
    # margin       How much extra space to use.
    import graphconst as gc
    
    fx, fy, fz, fw, fh, fd = object
    lw, lh, ld = label
    if alternate:
        lw, lh, ld = altlabel
    if orientation == gc.LABEL_VERTICAL:
        lw, lh = lh, lw
    if position == gc.LABEL_TOP:
        x = (fx+fw/2.0)-lw/2.0
        y = fy - lh - margin
    elif position == gc.LABEL_RIGHT:
        x = fx + fw + margin
        y = (fy+fh/2.0)-lh/2.0
    elif position == gc.LABEL_BOTTOM:
        x = (fx+fw/2.0)-lw/2.0
        y = fy + fh + margin
    elif position == gc.LABEL_LEFT:
        x = fx - lw - margin
        y = (fy+fh/2.0)-lh/2.0
    elif position == gc.LABEL_TOP_LEFT:
        x = fx - lw - margin
        y = fy - lh - margin
    elif position == gc.LABEL_TOP_RIGHT:
        x = fx + fw + margin
        y = fy - lh - margin
    elif position == gc.LABEL_BOTTOM_LEFT:
        x = fx - lw - margin
        y = fy + fh + margin
    elif position == gc.LABEL_BOTTOM_RIGHT:
        x = fx + fw + margin
        y = fy + fh + margin
    else:
        raise AssertionError, "Unknown position"
    z = fz + fd/2.0 - ld/2.0
    return x, y, z

def _score_label(
    index, label, altlabel, position, orientation, alternate,
    objects, labels, altlabels, coords, positions, orientations, alternates, 
    bounds, margin, overlap_cache):
    import graphconst as gc
    assert len(objects) == len(labels)
    assert len(labels) == len(coords)
    w, h, d = label
    if alternate:
        w, h, d = altlabel
    if orientation == gc.LABEL_VERTICAL:
        w, h = h, w
    
    coord = _calc_label_coord(
        objects[index], label, altlabel, position, orientation, alternate,
        margin)
    x, y, z = coord

    # Use a scoring heuristic that takes into account:
    # 1.  The position of the label.  Some are favored over others.
    # 2.  If the label overlaps with other labels.
    # 3.  If the label overlaps with objects on the plot.
    # 4.  If the label is out of bounds.
    #
    # Tried both all-or-nothing penalties and penalties that take into
    # account the amount of overlap.  The ones that give partial
    # credit generates much better solutions.
    PENALTY_ALTERNATE = -10   # Whether to use the alternate label.
    PENALTY_LABEL = -400      # Overlaps with a label.
    PENALTY_OBJECT = -200     # Overlaps with an object.
    PENALTY_OUTSIDE = -10000  # Outside the boundary.
    
    score = 0

    # Penalize if using the alternate label.
    if alternate:
        score += PENALTY_ALTERNATE

    # Score based on the position.
    position2score = {
        gc.LABEL_TOP : 0,
        gc.LABEL_BOTTOM : 0,
        gc.LABEL_LEFT : -5,
        gc.LABEL_RIGHT : -5,
        gc.LABEL_TOP_LEFT : -20,
        gc.LABEL_TOP_RIGHT : -20,
        gc.LABEL_BOTTOM_LEFT : -20,
        gc.LABEL_BOTTOM_RIGHT : -20,
        }

    favor_horizontal = True
    fx, fy, fz, fw, fh, fd = objects[index]
    if fh > fd*2:
        # If the object to be labelled is long and skinny, then a
        # vertical label would look better.
        favor_horizontal = False

    pos_score = position2score[position]
    if favor_horizontal and orientation == gc.LABEL_VERTICAL:
        pos_score -= 10   # favor horizontal ones.
    if not favor_horizontal and orientation == gc.LABEL_HORIZONTAL:
        pos_score -= 10   # favor vertical ones.
    score += pos_score

    # Penalize overlaps with objects.
    key = index, position, orientation, alternate
    I = overlap_cache.get(key, range(len(objects)))
    #I = range(len(objects))
    for i in I:
        fx, fy, fz, fw, fh, fd = objects[i]
        p = _perc_overlap((x, y, z, w, h, d), objects[i])
        #if p:
        #    score += PENALTY_OBJECT
        score += p*PENALTY_OBJECT
        #if i in [114] and orientation == 0 and position == 2:
        #    xx = ["HERE2", i, position, orientation, score] + \
        #         list(objects[i])[:2] + list((x, y, w, h))
        #    print "\t".join(map(str, xx))

    # Penalize overlaps with other labels.
    for i in range(len(labels)):
        if i == index:
            continue
        lx, ly, lz = coords[i]
        lw, lh, ld = labels[i]
        if alternates[i]:
            lw, lh, ld = altlabels[i]
        if (lw, lh, ld) == (0, 0, 0):
            continue
        if orientations[i] == gc.LABEL_VERTICAL:
            lw, lh = lh, lw
        p = _perc_overlap((x, y, z, w, h, d), (lx, ly, lz, lw, lh, ld))
        #if p:
        #    score += PENALTY_LABEL
        score += p*PENALTY_LABEL
            
    # Check if this label is in bounds.
    bx, by, bz, bw, bh, bd = bounds
    ib_x = x >= bx and x+w < bx+bw
    ib_y = y >= by and y+h < by+bh
    #ib_z = z >= bz and z+d < bz+bd
    ib_z = True  # don't care about z coordinate.
    if not (ib_x and ib_y and ib_z):
        score += PENALTY_OUTSIDE
    return score

def _position_one_label(
    index, label, altlabel, 
    objects, labels, altlabels, coords, positions, orientations, alternates, 
    bounds, margin, overlap_cache):
    import itertools
    import graphconst as gc

    assert len(objects) == len(labels)
    assert len(labels) == len(coords)

    # Find the best position for the label.
    ALL_POSITIONS = [
        gc.LABEL_TOP, gc.LABEL_BOTTOM, gc.LABEL_LEFT, gc.LABEL_RIGHT,
        gc.LABEL_TOP_LEFT, gc.LABEL_TOP_RIGHT,
        gc.LABEL_BOTTOM_LEFT, gc.LABEL_BOTTOM_RIGHT]
    ALL_ORIENTATIONS = [
        gc.LABEL_HORIZONTAL, gc.LABEL_VERTICAL]
    ALL_ALTERNATES = [False, True]
    
    best_pos = best_orient = best_alt = max_score = None
    for x in itertools.product(
        ALL_POSITIONS, ALL_ORIENTATIONS, ALL_ALTERNATES):
        pos, orient, alt = x
        score = _score_label(
            index, label, altlabel, pos, orient, alt, 
            objects, labels, altlabels, coords, positions, orientations,
            alternates, bounds, margin, overlap_cache)
          
        #x = "HERE1", index, pos, orient, alt, score
        #print "\t".join(map(str, x))
        if max_score is None or score > max_score:
            best_pos, best_orient, best_alt, max_score = (
                pos, orient, alt, score)
    return best_pos, best_orient, best_alt

def _cache_label_object_overlaps(objects, labels, altlabels, margin):
    import itertools
    import graphconst as gc

    # Optimization: cache the overlaps between labels and objects.
    ALL_POSITIONS = [
        gc.LABEL_TOP, gc.LABEL_BOTTOM, gc.LABEL_LEFT, gc.LABEL_RIGHT,
        gc.LABEL_TOP_LEFT, gc.LABEL_TOP_RIGHT,
        gc.LABEL_BOTTOM_LEFT, gc.LABEL_BOTTOM_RIGHT]
    ALL_ORIENTATIONS = [
        gc.LABEL_HORIZONTAL, gc.LABEL_VERTICAL]
    ALL_ALTERNATES = [False, True]

    overlaps = {}  # i, pos, orient, alternate -> list of objects
    for x in itertools.product(
        range(len(labels)), ALL_POSITIONS, ALL_ORIENTATIONS, ALL_ALTERNATES):
        i, pos, orient, alt = x
        if alt and altlabels[i] == (0, 0, 0):
            continue
        if not alt and labels[i] == (0, 0, 0):
            continue
        overlaps[(i, pos, orient, alt)] = []
        coord = _calc_label_coord(
            objects[i], labels[i], altlabels[i], pos, orient, alt, margin)
        x, y, z = coord
        w, h, d = labels[i]
        if alt:
            w, h, d = altlabels[i]
        if orient == gc.LABEL_VERTICAL:
            w, h = h, w
        for j in range(len(objects)):
            if _overlaps((x, y, z, w, h, d), objects[j]):
                overlaps[(i, pos, orient, alt)].append(j)
    return overlaps

def position_labels(bounds, objects, labels, altlabels, margin):
    # bounds is (x, y, z, width, height, depth) for the whole plotting
    # area.  objects is a list of (x, y, z, width, height, depth) of
    # the objects on the plot to avoid.  labels is a list of the
    # (width, height, depth) for the label of each object.  labels
    # should be parallel to objects.  If an object is not labeled,
    # then width, height, and depth should each be 0.  altlabels is a
    # list of alternate labels.  It should be in the same form as
    # labels.  margin is the amount of space to add between the edge
    # of the object and the label.
    import graphconst as gc

    assert len(bounds) == 6
    assert len(objects) == len(labels)
    assert len(objects) == len(altlabels)
    for x in objects:
        assert len(x) == 6
    for x in labels:
        assert len(x) == 3
    for x in altlabels:
        assert len(x) == 3

    # Initialize the positions and coordinates for each label.
    positions = [gc.LABEL_TOP] * len(labels)
    orientations = [gc.LABEL_HORIZONTAL] * len(labels)
    alternates = [False] * len(labels)
    coords = [None] * len(labels)
    for i in range(len(labels)):
        coords[i] = _calc_label_coord(
            objects[i], labels[i], altlabels[i], positions[i], orientations[i],
            alternates[i], margin)

    # Cache overlapping objects and the labels to prevent a long
    # search each time.
    overlap_cache = _cache_label_object_overlaps(
        objects, labels, altlabels, margin)

    MAX_ITER = 50
    MAX_SAME = 2   # how many times same labels moved back and forth.
    
    num_iter = 0
    num_same_moved = 0
    moved = None
    while num_same_moved < MAX_SAME and num_iter < MAX_ITER:
        num_iter += 1

        # For each label, find the best position.
        last_moved = moved
        moved = []
        for i in range(len(labels)):
            if labels[i] == (0, 0, 0):
                continue
            # Check the positioning of this label.  If it's good, then
            # don't move it.
            score = _score_label(
                i, labels[i], altlabels[i], positions[i], orientations[i],
                alternates[i], objects, labels, altlabels, coords,
                positions, orientations, alternates, bounds, margin,
                overlap_cache)
            if score >= 0:
                continue

            # Try to find a better position for this label.
            x = _position_one_label(
                i, labels[i], altlabels[i],
                objects, labels, altlabels,
                coords, positions, orientations, alternates,
                bounds, margin, overlap_cache)
            pos, orient, alt = x
            if(pos == positions[i] and orient == orientations[i] and
               alt == alternates[i]):
                # Hasn't moved, don't do anything.
                continue
            
            positions[i] = pos
            orientations[i] = orient
            alternates[i] = alt
            coords[i] = _calc_label_coord(
                objects[i], labels[i], altlabels[i],
                positions[i], orientations[i], alternates[i], margin)
            moved.append(i)
        if last_moved == moved:
            num_same_moved += 1
        else:
            num_same_moved = 0
        #print moved, num_same_moved, num_iter; import sys; sys.stdout.flush()
    return positions, orientations, alternates, coords


def test_find_heatmap_size():
    print find_tall_heatmap_size(10, 10)  # 12, 20
    print find_tall_heatmap_size(5, 10)   # 6, 20
    print find_tall_heatmap_size(10, 10, min_box_width=100)  # 100, 161
    print find_tall_heatmap_size(5, 10, min_box_width=100)   # 100, 323

    print find_wide_heatmap_size(10, 10)  # 20, 12
    print find_wide_heatmap_size(5, 10)   # 16, 20

def test_place_ticks():
    print place_ticks(0.00, 1.00, 5)
    print place_ticks(0.03, 0.94, 3)
    print place_ticks(-12, 15, 10)

def test_choose_tick_delta():
    print _choose_tick_delta(-2.0, 2.0, 6)
    print _choose_tick_delta(0, 1)
    print _choose_tick_delta(0, 100)
    print _choose_tick_delta(0, 21)
    print _choose_tick_delta(0, 0.8, 2)
    print _choose_tick_delta(0, 1, 40)
    print _choose_tick_delta(-12, 15, 2)
    print _choose_tick_delta(0, 99, 2)
    print _choose_tick_delta(0.0, 0.1, 10)

if __name__ == '__main__':
    #test_find_heatmap_size()
    test_choose_tick_delta()
    #test_place_ticks()
