"""Cross-platform Tk window icon for the ground-truthing GUI."""

from __future__ import annotations

from pathlib import Path
import tkinter as tk

from nps_active_space import ACTIVE_SPACE_DIR

_ICON_ATTR = "_nps_window_icon_photo"


def nps_logo_paths() -> tuple[Path, Path]:
    img_dir = Path(ACTIVE_SPACE_DIR) / "img"
    return img_dir / "flat-four-color.png", img_dir / "flat-four-color.ico"


def apply_window_icon(window: tk.Misc) -> None:
    """Set the NPS flat-four logo when packaged assets are present.

    Linux Tk often rejects ``iconbitmap`` for ``.ico`` files; ``iconphoto`` with PNG
    works on X11 and Windows. Missing assets are ignored so the GUI still launches.
    """
    png_path, ico_path = nps_logo_paths()
    if png_path.is_file():
        try:
            photo = tk.PhotoImage(file=str(png_path))
            window.iconphoto(True, photo)
            setattr(window, _ICON_ATTR, photo)
            return
        except tk.TclError:
            pass
    if ico_path.is_file():
        try:
            window.iconbitmap(str(ico_path))
        except tk.TclError:
            pass
