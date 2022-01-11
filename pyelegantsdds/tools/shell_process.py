import shlex
import subprocess as subp


def call(cmdstr, rootname):
    # general script to execute a shell command, and read the stderr from
    # the underlying executable root name
    
    errors, output = None, None
    try:
        task = subp.run(shlex.split(cmdstr), check=True, shell=False, stdout=subp.PIPE,
                           stderr=subp.STDOUT)
        errors = task.stderr
        output = task.stdout
    except subp.CalledProcessError as e:
        print ('>>> AN ERROR HAS OCCURED! <<<')
        print (f'Check {rootname}.stderr')
        errors = e.output

    if output == None:
        output = ''
    else:
        output = output.decode('utf-8')
    with open(f'{rootname}.stdout', 'w') as f:
        f.write(output)

    if errors == None:
        errors = ''
    else:
        errors = errors.decode('utf-8')
    with open(f'{rootname}.stderr', 'w') as f:
        f.write(errors)