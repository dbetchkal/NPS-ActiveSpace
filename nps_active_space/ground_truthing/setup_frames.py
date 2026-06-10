import tkinter as tk
from PIL import Image, ImageTk

from nps_active_space import ACTIVE_SPACE_DIR
from nps_active_space.ground_truthing.base import _AppFrame
from nps_active_space.ground_truthing.session_frames import _AnnotationLoadFrame


class _WelcomeFrame(_AppFrame):
    """
    The opening frame for the Ground Truthing application that welcomes the user.

    Parameters
    ----------
    master : tk.Tk
        The tkinter window this frame will be shown in.
    """
    def __init__(self, master):
        super().__init__(master)

        # Define widgets.
        frame_label = tk.Label(
            self,
            text='Welcome to the NPS Active Space Project\nGround Truthing Module!',
            font=('Avenir', 20, 'bold'),
            bg='ivory2'
        )
        continue_button = tk.Button(
            self,
            text='Continue >>',
            width=20,
            font=('Avenir', 8),
            bg='ivory2',
            command=lambda: self.master.switch_frame(_AnnotationLoadFrame)
        )
        im = Image.open(f"{ACTIVE_SPACE_DIR}/img/flat-four-color.png").resize((138, 181))
        nps_logo = ImageTk.PhotoImage(im)
        label = tk.Label(self, image=nps_logo, bg='ivory2')
        label.image = nps_logo  # NOTE: This re-definition is required for windows machines.

        # Place widgets.
        label.place(relx=0.5, rely=0.3, anchor='center')
        frame_label.place(relx=0.5, rely=0.55, anchor='center')
        continue_button.place(relx=0.9, rely=0.9, anchor='center')
