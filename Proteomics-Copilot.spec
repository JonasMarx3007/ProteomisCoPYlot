# -*- mode: python ; coding: utf-8 -*-

import os
from PyInstaller.utils.hooks import copy_metadata
import streamlit

project_dir = os.path.abspath(os.getcwd())

datas = [
    (os.path.join(project_dir, "app.py"), "."),
    (os.path.join(project_dir, "ui"), "ui"),
    (os.path.join(project_dir, "utils"), "utils"),
    (os.path.join(project_dir, "favicon.ico"), "."),
] + copy_metadata('streamlit')

streamlit_path = os.path.dirname(streamlit.__file__)
static_folder = os.path.join(streamlit_path, 'static')

for root, dirs, files in os.walk(static_folder):
    for file in files:
        full_path = os.path.join(root, file)
        rel_path = os.path.relpath(root, streamlit_path)
        datas.append((full_path, os.path.join('streamlit', rel_path)))

a = Analysis(
    ['run_app.py'],
    pathex=[project_dir],
    binaries=[],
    datas=datas,
    hiddenimports = [
    'streamlit.runtime.scriptrunner.magic_funcs',
    'matplotlib',
    'matplotlib.backends.backend_tkagg',
    'matplotlib.backends.backend_agg',
    'scipy',
    'scipy.spatial',
    'scipy.sparse',
    'scipy.linalg',
    'scipy.optimize',
    'scipy.stats',
    'sklearn',
    'sklearn.decomposition',
    'sklearn.utils._cython_blas',
    'statsmodels.stats.multitest',
    ],
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=[],
    noarchive=False,
    optimize=0,
)

pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='Proteomics-Copilot',
    debug=False,                        #Debug
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,                       # Console
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon='receptor_icon.ico',
)
