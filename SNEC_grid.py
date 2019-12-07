import asyncio
import subprocess

async def snec(dir):
    subprocess.run(['cp', '/home/sbracha/SNEC/SNEC-1.01/'+ dir + '/parameters', './'], universal_newlines=True)
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





asyncio.run(snec('Ni014_E210_Mix6'))
asyncio.run(snec('Ni014_E210_Mix9'))
asyncio.run(snec('Ni014_E210_Mix12'))

asyncio.run(snec('Ni014_E300_Mix6'))
asyncio.run(snec('Ni014_E300_Mix9'))
asyncio.run(snec('Ni014_E300_Mix12'))

asyncio.run(snec('Ni016_E210_Mix6'))
asyncio.run(snec('Ni016_E210_Mix9'))
asyncio.run(snec('Ni016_E210_Mix12'))

asyncio.run(snec('Ni016_E300_Mix6'))
asyncio.run(snec('Ni016_E300_Mix9'))
asyncio.run(snec('Ni016_E300_Mix12'))




subprocess.run(['cp', 'Ni014_E210_Mix6/parameters', './'], universal_newlines=True)
subprocess.run(['./snec'], universal_newlines=True)

subprocess.run(['cp', 'Ni014_E210_Mix9/parameters', './'], universal_newlines=True)
subprocess.run(['./snec'], universal_newlines=True)

subprocess.run(['cp', 'Ni014_E210_Mix12/parameters', './'], universal_newlines=True)
subprocess.run(['./snec'], universal_newlines=True)



subprocess.run(['cp', 'Ni014_E300_Mix6/parameters', './'], universal_newlines=True)
subprocess.run(['./snec'], universal_newlines=True)

subprocess.run(['cp', 'Ni014_E300_Mix9/parameters', './'], universal_newlines=True)
subprocess.run(['./snec'], universal_newlines=True)

subprocess.run(['cp', 'Ni014_E300_Mix12/parameters', './'], universal_newlines=True)
subprocess.run(['./snec'], universal_newlines=True)




subprocess.run(['cp', 'Ni016_E210_Mix6/parameters', './'], universal_newlines=True)
subprocess.run(['./snec'], universal_newlines=True)

subprocess.run(['cp', 'Ni016_E210_Mix9/parameters', './'], universal_newlines=True)
subprocess.run(['./snec'], universal_newlines=True)

subprocess.run(['cp', 'Ni016_E210_Mix12/parameters', './'], universal_newlines=True)
subprocess.run(['./snec'], universal_newlines=True)



subprocess.run(['cp', 'Ni016_E300_Mix6/parameters', './'], universal_newlines=True)
subprocess.run(['./snec'], universal_newlines=True)

subprocess.run(['cp', 'Ni016_E300_Mix9/parameters', './'], universal_newlines=True)
subprocess.run(['./snec'], universal_newlines=True)

subprocess.run(['cp', 'Ni016_E300_Mix12/parameters', './'], universal_newlines=True)
subprocess.run(['./snec'], universal_newlines=True)
