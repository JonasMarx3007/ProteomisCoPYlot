import os
import sys
import streamlit.web.cli as stcli
import contextlib
import io
import webbrowser
import time
import threading


def resource_path(relative_path):
    try:
        base_path = sys._MEIPASS
    except AttributeError:
        base_path = os.path.abspath(".")
    return os.path.join(base_path, relative_path)

def open_browser(url, delay=3):
    time.sleep(delay)
    webbrowser.open(url)

def main():
    os.environ["STREAMLIT_GLOBAL_DEVELOPMENT_MODE"] = "false"

    app_path = resource_path("app.py")
    sys.argv = [
        "streamlit", "run", app_path,
        "--server.port=8501",
        "--server.headless=true"
    ]

    BLUE = '\033[94m'
    RESET = '\033[0m'

    print("\n=== Proteomics-Copilot Streamlit App ===")
    print(f"The app is running at: {BLUE}http://localhost:8501{RESET}")
    print()
    print("To open the app in your browser:")
    print(" - On Windows: Hold Ctrl (or Cmd on Mac) and click the link, or copy-paste it manually.")
    print(" - On other systems: Copy-paste the link into your browser.")
    print()
    print("To stop the app, return here and press Ctrl+C (Cmd+C on Mac).")
    print("=========================================\n")

    threading.Thread(target=open_browser, args=("http://localhost:8501",), daemon=True).start()

    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        sys.exit(stcli.main())

if __name__ == "__main__":
    main()
