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

    if (dataset_dir / 'molprobity.out').exists():
        print(f'Already have molprobity results at {dataset_dir}!')
        return
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
    dataset_dirs = []
    for system_dir in pathlib.Path(path).glob('*'):
        if not system_dir.is_dir():
            continue
        print(system_dir.name)
        for dataset_dir in (system_dir / 'processed_datasets').glob('*'):
            if not dataset_dir.is_dir():
                continue
            if  (dataset_dir / f'{dataset_dir.name}-pandda-input.pdb').exists():
                dataset_dirs.append(dataset_dir)
    print(f'Got {len(dataset_dirs)}')
    run_molprobity(dataset_dirs[0])

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