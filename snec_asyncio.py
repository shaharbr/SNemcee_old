import asyncio, os

async def snec(dir):
    os.chdir(os.path.abspath('/home/sbracha/SNEC/SNEC_'+dir))
    cmd = '/home/sbracha/SNEC/SNEC-1.01/snec'
    proc = await asyncio.create_subprocess_shell(
        cmd,
        stdout=asyncio.subprocess.PIPE,
        stderr=asyncio.subprocess.PIPE)
    stdout, stderr = await proc.communicate()
    print(f'[{cmd!r} exited with {proc.returncode}]')
    if stdout:
        print(f'[stdout]\n{stdout.decode()}')
    if stderr:
        print(f'[stderr]\n{stderr.decode()}')


async def main():
    await asyncio.gather(
        snec('aNi014_E240_Mix1_Mzams21'),
        )



asyncio.run(main())
