import subprocess
import pathlib
import os
import time
import glob
import datetime

project_dir = pathlib.Path(os.environ['CONDA_PREFIX']).parent
database = pathlib.Path(project_dir, 'database')
exe_dir = pathlib.Path(project_dir, 'bin')
src_dir = pathlib.Path(project_dir, 'sources')
subsrc_dir = [path for path in pathlib.Path(project_dir, 'sources').iterdir() if path.is_dir()]

exe_name = pathlib.Path(__file__).stem
path_main = pathlib.Path(exe_dir, f'{str(exe_name)}.cpp')
exe = pathlib.Path(exe_dir, path_main.stem)

input_dir = pathlib.Path('/data/HiRA_Cali')
out_dir = pathlib.Path('/data/kin/hira')


def main():

    # compiling the program
    root_inc = subprocess.run('root-config --cflags --libs --glibs',
                              capture_output=True, shell=True, text=True, encoding='utf-8')
    root_inc = root_inc.stdout.strip()
    
    src_inc = ' '.join([f'-I{str(src_dir)}'] + [f'-I{str(path)}' for path in subsrc_dir] + [f'-I{root_inc}'])

    print('compiling...\n')
    compile_cmd = f'g++ {str(path_main)} -o {str(exe)} {src_inc}'
    print(compile_cmd, '\n')
    subprocess.run(compile_cmd, shell=True, text=True)

    runs('Ca', 'Ni', 40, 58, 56)
    runs('Ca', 'Ni', 40, 58, 140)
    runs('Ca', 'Ni', 48, 64, 56)
    runs('Ca', 'Ni', 48, 64, 140)

    runs('Ca', 'Sn', 40, 112, 56)
    runs('Ca', 'Sn', 40, 112, 140)
    runs('Ca', 'Sn', 48, 124, 56)
    runs('Ca', 'Sn', 48, 124, 140)



def runs(beam, target, beamA, targetA, energy):
    indices = get_indices(beam, target, beamA, targetA, energy)
    for id in indices:
        run(beam, target, beamA, targetA, energy, id)


def get_indices(beam, target, beamA, targetA, energy):
    path_data = pathlib.Path(
        input_dir, f'{beamA}{beam}{targetA}{target}_{energy}MeVu')
    
    files = glob.glob(f'{str(path_data)}/CalibratedData_*.root')
    indices = [pathlib.Path(f).stem[-4:] for f in files]

    return indices

def run(beam, target, beamA, targetA, energy, id):
    path_data = pathlib.Path(
        input_dir, f'{beamA}{beam}{targetA}{target}_{energy}MeVu/CalibratedData_{id}.root')

    sub_dir = pathlib.Path(out_dir, f'{beamA}{beam}{targetA}{target}_{energy}MeVu')
    sub_dir.mkdir(exist_ok=True, parents=True)
    path_out = pathlib.Path(sub_dir, f'run-{id}.root')

    args = f'{str(path_data)} {str(path_out)}'
    start_time = time.time()
    subprocess.run(f'{str(exe)} {args}', shell=True, text=True)
    finish_time = time.time()
    elapsed_time = finish_time - start_time
    print(f'reaction {beam}{beamA}{target}{targetA}E{energy}-run{id} done in {datetime.timedelta(elapsed_time)}')



if __name__ == '__main__':
    main()
