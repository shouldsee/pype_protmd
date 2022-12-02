from pype import Controller
from pype import check_write_single_target as CWST

def know_puff(ctl):
    ctl.lazy_wget('http://boscoh.com/puff/puff-v1.0.zip',name='download_puff')
    x = 'puff-v1.0.zip'
    ctl.RWC( CWST, x+'.done', run=f'''unzip {x}; echo 1 >{x}.done''')
#    ctl.RWC(run='2to3 ../puff/package/pdbtool/ -wno ../puff/package/pdbtool.py3')
    ctl.RWC(CWST, 'package/pdbtool.py3/puff.py', 
    run='''
    rm -rf package/pdbtool_py3; 
    cp -rf package/pdbtool package/pdbtool_py3;
    python3 -m lib2to3 -Wn ../puff/package/pdbtool_py3
    ''')

    '''
    rather extensive modifications
    [todo] modify puff to include
    '''
    ### manual changes
    ### vector3d indentation error

#    ctl.RWC(run='python3 my223.py -f py223_fixers ../puff/package/pdbtool/ -Wn -o ../puff/package/pdbtool.py3')
    #ctl.cd('package/pdbtool/')
#    ctl.RWC(run='python3 -mpdb -cc package/pdbtool.py3/puff.py package/example')
    ctl.RWC(run='rm -rf package/example/hairpin/* ')
    ctl.RWC(run='python3 package/pdbtool.py3/puff.py package/example && echo done')
#from pype import Controller
if __name__ == '__main__':
    x = Controller()
    know_puff(x)
    x.build('$HOME/catsmile/prot/puff')
    pass

