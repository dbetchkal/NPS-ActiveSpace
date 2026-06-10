from matplotlib.widgets import RangeSlider


class FastRangeSlider(RangeSlider):
    """
    Subclass of RangeSlider for compatibility with blitting.
    RangeSlider is a widget, not an Artist, so it doesn't have a draw method.
    This causes issues when we want to use ax.draw_artist() for blitting compatibility.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Set all changing UI components to be animated, so they don't render by default.
        # This allows us to capture the background without these UI components already there.
        self.poly.set_animated(True)
        for h in self._handles:
            h.set_animated(True)

    def draw(self):
        # call draw_artist() on each changing UI component individually
        self.ax.draw_artist(self.poly)
        for h in self._handles:
            self.ax.draw_artist(h)
