import tkinter as tk
from abc import ABC


class _AppFrame(ABC, tk.Frame):
    """
    Abstract base class for all application frames.

    Parameters
    ----------
    master : tk.Tk
        A tkinter app instance that will display the frame.
    """
    def __init__(self, master: tk.Tk) -> None:
        super().__init__(master)
        self.master = master
