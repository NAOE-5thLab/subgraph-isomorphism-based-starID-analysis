for /l %%a in (2,1,10) do (
    for /l %%b in (0,1,10) do (
        for /l %%c in (0,1,4) do (
            start python subgraph_monte.py --n_obs %%a --theta_img %%b --k %%c --Vmax 0 --theta_FOV 0
            start python subgraph_monte.py --n_obs %%a --theta_img %%b --k %%c --Vmax 0 --theta_FOV 1
            start python subgraph_monte.py --n_obs %%a --theta_img %%b --k %%c --Vmax 0 --theta_FOV 2
            start python subgraph_monte.py --n_obs %%a --theta_img %%b --k %%c --Vmax 0 --theta_FOV 3
            start python subgraph_monte.py --n_obs %%a --theta_img %%b --k %%c --Vmax 1 --theta_FOV 0
            start python subgraph_monte.py --n_obs %%a --theta_img %%b --k %%c --Vmax 1 --theta_FOV 1
            start python subgraph_monte.py --n_obs %%a --theta_img %%b --k %%c --Vmax 1 --theta_FOV 2
            start python subgraph_monte.py --n_obs %%a --theta_img %%b --k %%c --Vmax 1 --theta_FOV 3
            start python subgraph_monte.py --n_obs %%a --theta_img %%b --k %%c --Vmax 2 --theta_FOV 0
            start python subgraph_monte.py --n_obs %%a --theta_img %%b --k %%c --Vmax 2 --theta_FOV 1
            start python subgraph_monte.py --n_obs %%a --theta_img %%b --k %%c --Vmax 2 --theta_FOV 2
            start python subgraph_monte.py --n_obs %%a --theta_img %%b --k %%c --Vmax 2 --theta_FOV 3
            start python subgraph_monte.py --n_obs %%a --theta_img %%b --k %%c --Vmax 3 --theta_FOV 0
            start python subgraph_monte.py --n_obs %%a --theta_img %%b --k %%c --Vmax 3 --theta_FOV 1
            start python subgraph_monte.py --n_obs %%a --theta_img %%b --k %%c --Vmax 3 --theta_FOV 2
            start python subgraph_monte.py --n_obs %%a --theta_img %%b --k %%c --Vmax 3 --theta_FOV 3
            start python subgraph_monte.py --n_obs %%a --theta_img %%b --k %%c --Vmax 4 --theta_FOV 0
            start python subgraph_monte.py --n_obs %%a --theta_img %%b --k %%c --Vmax 4 --theta_FOV 1
            start python subgraph_monte.py --n_obs %%a --theta_img %%b --k %%c --Vmax 4 --theta_FOV 2
            start /wait python subgraph_monte.py --n_obs %%a --theta_img %%b --k %%c --Vmax 4 --theta_FOV 3
        )
    )
)

