import subprocess
import pathlib
import os
import time

project_dir = pathlib.Path(os.environ['CONDA_PREFIX']).parent
database = pathlib.Path(project_dir, 'database')
exe_dir = pathlib.Path(project_dir, 'analysis')
src_dir = pathlib.Path(project_dir, 'sources')
subsrc_dir = [path for path in pathlib.Path(project_dir, 'sources').iterdir() if path.is_dir()]

exe_name = pathlib.Path(__file__).stem
path_main = pathlib.Path(exe_dir, f'{str(exe_name)}.cpp')

input_dir = pathlib.Path('/data/HiRA_Cali')
out_dir = pathlib.Path(project_dir, 'result/spectra')
out_dir.mkdir(exist_ok=True, parents=True)

beam = 'Ca'
target = 'Ni'
beamA = 48
targetA = 64
energy = 140
bimp = (0., 0.4)


def main():
    root_inc = subprocess.run('root-config --cflags --libs --glibs',
                              capture_output=True, shell=True, text=True, encoding='utf-8')
    root_inc = root_inc.stdout.strip()
    exe = pathlib.Path(exe_dir, path_main.stem)

    src_inc = ' '.join([f'-I{str(src_dir)}'] + [f'-I{str(path)}' for path in subsrc_dir] + [f'-I{root_inc}'])

    print('compiling...\n')
    compile_cmd = f'g++ {str(path_main)} -o {str(exe)} {src_inc}'
    print(compile_cmd, '\n')
    subprocess.run(compile_cmd, shell=True, text=True)

    path_data = pathlib.Path(
        input_dir, f'{beamA}{beam}{targetA}{target}_{energy}MeVu')

    path_out = pathlib.Path(
        out_dir, f'{beam}{beamA}{target}{targetA}E{energy}_bmin{bimp[0]:.1f}_bmax{bimp[1]:.1f}.root')


    path_runinfo = pathlib.Path(database, 'config/RunInfo.data')
    path_badmap = pathlib.Path(database, 'GeoEff')
    path_angles = pathlib.Path(database, 'Cal_PixelAngle/PixelAngle_BeamPos_0_0_0.dat')

    inputs = list(
        map(str, [f'{beam}{beamA}{target}{targetA}E{energy}', path_data, path_out, *bimp, path_runinfo, path_badmap, path_angles]))
    args = ' '.join(inputs)
    print('Running...\n')
    print(f'{str(exe)} {args}')

    start_time = time.time()
    subprocess.run(f'{str(exe)} {args}', shell=True, text=True)
    finish_time = time.time()
    elapsed_time = finish_time - start_time
    print(f'All processes done in {elapsed_time}')




if __name__ == '__main__':
    main()
