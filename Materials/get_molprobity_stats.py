import joblib
import subprocess
import fire
import pathlib

SCRIPT = '''module load ccp4\n
cd {dataset_dir}\n
molprobity.molprobity {pdb_path} {mtz_path}
'''

def run_molprobity(dataset_dir):
    dtag = dataset_dir.name
    st_path = dataset_dir / f'{dtag}-pandda-input.pdb'
    mtz_path = dataset_dir / f'{dtag}-pandda-input.mtz'
    script = SCRIPT.format(
        dataset_dir=dataset_dir,
        pdb_path=st_path,
        mtz_path = mtz_path
    )
    print(script)
    p=subprocess.Popen(
        script,
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    stdout, stderr = p.communicate()
    print(stdout)
    print(stderr)
    ...

def main(path):
    dataset_dirs = [dataset_dir for system_dir in pathlib.Path(path).glob('*') for dataset_dir in system_dir]
    print(f'Got {len(dataset_dirs)}')
    exit()
    
    joblib.Parallel(n_jobs=-1, verbose=50)(
        joblib.delayed(
            run_molprobity
            )(
                dataset_dir
            )
            for dataset_dir
            in dataset_dirs
        )
    
    ...

if __name__ == "__main__":
    fire.Fire(main)