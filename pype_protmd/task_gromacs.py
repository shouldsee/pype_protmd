'''
REF: http://www.mdtutorials.com/gmx/lysozyme/02_topology.html
# TODO: GRACE fails silently???
Runtime: Approx 2 mins on a RTX3080
'''
from .depend_mol import know_gromacs
from pype import Controller,s
# from .depend_mol import check_write_single_target
from pype import FRAME,fstr
from pype import RuntimeObject as RO
from pype import AppendTypeChecker as ATC
from pype import check_write_single_target as CWST





'''
# http://localhost:6006/mdsrv/webapp/?struc=/mdsrv/file/1PGB/500000/md_0_1.gro&traj=file://1PGB/500000/md_0_1.xtc&mdsrv=/mdsrv/&data=/mdsrv/

Command line:
  gmx energy -f nvt.edr -o temperature.xvg

Opened nvt.edr as single precision energy file

Select the terms you want from the following list by
selecting either (part of) the name or the number or a combination.
End your selection with an empty line or a zero.
-------------------------------------------------------------------
  1  Bond             2  U-B              3  Proper-Dih.      4  Improper-Dih.
  5  CMAP-Dih.        6  LJ-14            7  Coulomb-14       8  LJ-(SR)
  9  Disper.-corr.   10  Coulomb-(SR)    11  Coul.-recip.    12  Position-Rest.
 13  Potential       14  Kinetic-En.     15  Total-Energy    16  Conserved-En.
 17  Temperature     18  Pres.-DC        19  Pressure        20  Constr.-rmsd
 21  Vir-XX          22  Vir-XY          23  Vir-XZ          24  Vir-YX
 25  Vir-YY          26  Vir-YZ          27  Vir-ZX          28  Vir-ZY
 29  Vir-ZZ          30  Pres-XX         31  Pres-XY         32  Pres-XZ
 33  Pres-YX         34  Pres-YY         35  Pres-YZ         36  Pres-ZX
 37  Pres-ZY         38  Pres-ZZ         39  #Surf*SurfTen   40  T-Protein
 41  T-non-Protein                       42  Lamb-Protein
 43  Lamb-non-Protein
'''

# from pype import FRAME, RuntimeObject, THIS


def build(ctl, GMX, PDB_ID, md_steps):
    '''
    performs type check before start building
    '''
    ctl.runtime_initer('GMX',  GMX, str)
    ctl.runtime_initer('PDB_ID',PDB_ID,str)
    ctl.runtime_initer('md_steps',md_steps)

    # ref_t = 301
    ref_t = 350
    # ctl.init_cd(f'{{PDB_ID}}/{md_steps}')
    # ctl.runtime_setter['PDB_ID'] = ATC(PDB_ID,str)

    def lazy_grace_png(TARGET):
        '''
        call grace to generate png
        '''
        SRC = TARGET[:-len('.png')]
        return ctl.RWC(CWST, check_ctx=TARGET, run=f'grace -nxy {SRC} -hdevice PNG -hardcopy -printfile {TARGET}')

    ctl.lazy_grace_png = lazy_grace_png

    def write_dep_file(rt):
        with open('ions.mdp','w') as f:
            f.write('''
; ions.mdp - used as input into grompp to generate ions.tpr
; Parameters describing what to do, when to stop and what to save
integrator  = steep         ; Algorithm (steep = steepest descent minimization)
emtol       = 1000.0        ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01          ; Minimization step size
nsteps      = 50000         ; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist         = 1         ; Frequency to update the neighbor list and long range forces
cutoff-scheme	= Verlet    ; Buffered neighbor searching
ns_type         = grid      ; Method to determine neighbor list (simple, grid)
coulombtype     = cutoff    ; Treatment of long range electrostatic interactions
rcoulomb        = 1.0       ; Short-range electrostatic cut-off
rvdw            = 1.0       ; Short-range Van der Waals cut-off
pbc             = xyz       ; Periodic Boundary Conditions in all 3 dimensions
            ''')
        with open('minim.mdp','w') as f:
            f.write('''
; minim.mdp - used as input into grompp to generate em.tpr
; Parameters describing what to do, when to stop and what to save
integrator  = steep         ; Algorithm (steep = steepest descent minimization)
emtol       = 1000.0        ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01          ; Minimization step size
nsteps      = 50000         ; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist         = 1         ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet    ; Buffered neighbor searching
ns_type         = grid      ; Method to determine neighbor list (simple, grid)
coulombtype     = PME       ; Treatment of long range electrostatic interactions
rcoulomb        = 1.0       ; Short-range electrostatic cut-off
rvdw            = 1.0       ; Short-range Van der Waals cut-off
pbc             = xyz       ; Periodic Boundary Conditions in all 3 dimensions
            ''')
        with open('nvt.mdp','w') as f:
            f.write('''
title                   = OPLS Lysozyme NVT equilibration
define                  = -DPOSRES  ; position restrain the protein
; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = 50000     ; 2 * 50000 = 100 ps
dt                      = 0.002     ; 2 fs
; Output control
nstxout                 = 500       ; save coordinates every 1.0 ps
nstvout                 = 500       ; save velocities every 1.0 ps
nstenergy               = 500       ; save energies every 1.0 ps
nstlog                  = 500       ; update log file every 1.0 ps
; Bond parameters
continuation            = no        ; first dynamics run
constraint_algorithm    = lincs     ; holonomic constraints
constraints             = h-bonds   ; bonds involving H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Nonbonded settings
cutoff-scheme           = Verlet    ; Buffered neighbor searching
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet
rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)
rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)
DispCorr                = EnerPres  ; account for cut-off vdW scheme
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling is on
tcoupl                  = V-rescale             ; modified Berendsen thermostat
tc-grps                 = Protein Non-Protein   ; two coupling groups - more accurate
tau_t                   = 0.1     0.1           ; time constant, in ps
ref_t                   = 300     300           ; reference temperature, one for each group, in K
; Pressure coupling is off
pcoupl                  = no        ; no pressure coupling in NVT
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Velocity generation
gen_vel                 = yes       ; assign velocities from Maxwell distribution
gen_temp                = 300       ; temperature for Maxwell distribution
gen_seed                = -1        ; generate a random seed
            ''')
        with open('npt.mdp','w') as f:
            f.write('''
title                   = OPLS Lysozyme NPT equilibration
define                  = -DPOSRES  ; position restrain the protein
; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = 50000     ; 2 * 50000 = 100 ps
dt                      = 0.002     ; 2 fs
; Output control
nstxout                 = 500       ; save coordinates every 1.0 ps
nstvout                 = 500       ; save velocities every 1.0 ps
nstenergy               = 500       ; save energies every 1.0 ps
nstlog                  = 500       ; update log file every 1.0 ps
; Bond parameters
continuation            = yes       ; Restarting after NVT
constraint_algorithm    = lincs     ; holonomic constraints
constraints             = h-bonds   ; bonds involving H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Nonbonded settings
cutoff-scheme           = Verlet    ; Buffered neighbor searching
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)
rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)
DispCorr                = EnerPres  ; account for cut-off vdW scheme
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling is on
tcoupl                  = V-rescale             ; modified Berendsen thermostat
tc-grps                 = Protein Non-Protein   ; two coupling groups - more accurate
tau_t                   = 0.1     0.1           ; time constant, in ps
ref_t                   = 300     300           ; reference temperature, one for each group, in K
; Pressure coupling is on
pcoupl                  = Parrinello-Rahman     ; Pressure coupling on in NPT
pcoupltype              = isotropic             ; uniform scaling of box vectors
tau_p                   = 2.0                   ; time constant, in ps
ref_p                   = 1.0                   ; reference pressure, in bar
compressibility         = 4.5e-5                ; isothermal compressibility of water, bar^-1
refcoord_scaling        = com
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Velocity generation
gen_vel                 = no        ; Velocity generation is off


            ''')
        with open('md.mdp','w') as f:
            # ref_t = 300
            f.write(f'''
title                   = OPLS Lysozyme production MD
; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = {md_steps}; 2 * 500000 = 1000 ps (1 ns)
;[DEBUG]
;nsteps                  = 10000    ; 2 * 5000 = 10 ps
dt                      = 0.002     ; 2 fs
; Output control
nstxout                 = 0         ; suppress bulky .trr file by specifying
nstvout                 = 0         ; 0 for output frequency of nstxout,
nstfout                 = 0         ; nstvout, and nstfout
nstenergy               = 5000      ; save energies every 10.0 ps
nstlog                  = 5000      ; update log file every 10.0 ps
nstxout-compressed      = 5000      ; save compressed coordinates every 10.0 ps
compressed-x-grps       = System    ; save the whole system
; Bond parameters
continuation            = yes       ; Restarting after NPT
constraint_algorithm    = lincs     ; holonomic constraints
constraints             = h-bonds   ; bonds involving H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Neighborsearching
cutoff-scheme           = Verlet    ; Buffered neighbor searching
ns_type                 = grid      ; search neighboring grid cells
nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)
rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling is on
tcoupl                  = V-rescale             ; modified Berendsen thermostat
tc-grps                 = Protein Non-Protein   ; two coupling groups - more accurate
tau_t                   = 0.1     0.1           ; time constant, in ps
ref_t                   = {ref_t} {ref_t}           ; reference temperature, one for each group, in K
; Pressure coupling is on
pcoupl                  = Parrinello-Rahman     ; Pressure coupling on in NPT
pcoupltype              = isotropic             ; uniform scaling of box vectors
tau_p                   = 2.0                   ; time constant, in ps
ref_p                   = 1.0                   ; reference pressure, in bar
compressibility         = 4.5e-5                ; isothermal compressibility of water, bar^-1
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Dispersion correction
DispCorr                = EnerPres  ; account for cut-off vdW scheme
; Velocity generation
gen_vel                 = no        ; Velocity generation is off
            ''')

    # CWST = check_write_single_target
    '''
    needs to be instantiate at build time
    '''

    #### build time dict is delayed until "build"
    #### buildtime is different from runtime in
    #### buildtime is state of
    # ctl.runtime_setter['GMX'] = 
    # ctl.runtime_setter['GMX'] = GMX = ctl.nodes["gromacs"].built+'/bin/gmx'
    # PDB_ID = ctl.runtime["PDB_ID"]


    ### getting variable from vain.
    '''
    PDB_ID needs to be instantiate at runtime
    '''
    # print(GMX())
    
    # ctl.lazy_wget( ctl.runtime.chain_with(lambda x:'''https://files.rcsb.org/download/{eprint(PDB_ID)}.pdb''',x)))
    ctl.lazy_wget( RO(ctl.runtime,'''https://files.rcsb.org/download/{PDB_ID}.pdb'''))
    # ctl.RWC(run=lambda x:ctl.runtime_mutable['PDB_ID']:=12))


    # ctl.RWC(run=lambda x:print(GMX.call()))
    ctl.RWC(run=write_dep_file)
    
    ### Clean versioned backup files
    ctl.RWC(run='rm -vf ./\#*')
    ### [more complex control flow to save a broken run]
    # ctl.RWC(run='rm -vf ./\#* ./*.gro ./*.top ')

    # ctl.RWC(run= lambda x:s(f'{GMX()} pdb2gmx -f {PDB_ID()}.pdb -o {PDB_ID()}_processed.gro -water spce -ff charmm27'))
    ctl.RWC(run= ('{GMX} pdb2gmx -f {PDB_ID}.pdb -o {PDB_ID}_processed.gro -water spce -ff charmm27'))

    ### adding periodic boundary
    '''
    -c centers
    -bt cubic boundary
    -d halved minimum distance between periodic images
    '''
    ctl.RWC(run='{GMX} editconf -f {PDB_ID}_processed.gro -o {PDB_ID}_newbox.gro -c -d 1.0')
    ctl.RWC(run='{GMX} solvate -cp {PDB_ID}_newbox.gro -o {PDB_ID}_solv.gro -cs spc216.gro -p topol.top')

    ctl.RWC(run='{GMX} grompp -f ions.mdp -c {PDB_ID}_solv.gro -p topol.top -o ions.tpr')
    ctl.RWC(run='echo 13 | {GMX} genion -s ions.tpr -o {PDB_ID}_solv_ions.gro -p topol.top -pname NA -nname CL -neutral')

    '''
    energy minimisation
    '''
    ctl.RWC(run='{GMX} grompp -o em.tpr -f minim.mdp -c {PDB_ID}_solv_ions.gro -p topol.top ')
    ctl.RWC(CWST,'em.gro',run='{GMX} mdrun -v -deffnm em')
    ctl.RWC(run='printf "10\n0\n" | {GMX} energy -f em.edr -o potential.xvg')
    ctl.lazy_grace_png( 'potential.xvg.png')

    '''
    equilibration
    '''
    ctl.RWC(run='''{GMX} grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr''')
    ctl.RWC(CWST,'nvt.gro',run='{GMX} mdrun -deffnm nvt')
    ctl.RWC(run='echo "17\n0\n" | {GMX} energy -f nvt.edr -o temperature.xvg')
    ctl.lazy_grace_png( 'temperature.xvg.png')

    ctl.RWC(run='''{GMX} grompp -f npt.mdp -c nvt.gro -r nvt.gro -p topol.top -o npt.tpr''')
    ctl.RWC(CWST,'npt.gro',run='{GMX} mdrun -deffnm npt')
#    ctl.RWC(ctx='npt.trr',run=f'{GMX} mdrun  -gputasks 0000  -nb gpu -pme gpu  -npme 1 -ntmpi 12 -deffnm npt')

    ctl.RWC(run='echo "19\n0\n" | {GMX} energy -f npt.edr -o pressure.xvg')
    ctl.RWC(run='echo "25\n0\n" | {GMX} energy -f npt.edr -o density.xvg')
    ctl.lazy_grace_png( 'density.xvg.png' )

    ctl.RWC(run='''{GMX} grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0_1.tpr''')
    ctl.RWC(CWST,'md_0_1.gro',run='{GMX} mdrun -deffnm md_0_1 -nb gpu')

    ctl.RWC(run=f'echo [SUCCESS]{__file__} 1>&2')
    return ctl



from pype import NotInitObject,ValueNotReadyError, print_tb_stacks
from pype import PlaceHolder

import argparse
def main():
    ctl = Controller()

    # ctl.runtime_setter['GMX'] = GMX = ctl.nodes["gromacs"].built+'/bin/gmx'
    # PDB_ID = PlaceHolder()
    pype1 = Controller.from_func(know_gromacs)
    # PDB_ID = pype1.runtime["PDB_ID"]
    PDB_ID = None
    '''
    Connecting pype 
    '''
    PDB_ID = PlaceHolder('PDB_ID')
    md_steps=5123457
        
    pype2 = Controller.from_func(build, 
        GMX    = pype1.built["GMX"],
        PDB_ID = PDB_ID.built,
        md_steps=md_steps,
        )
    # PDB_ID.set_pdb()
    parser= p = argparse.ArgumentParser()
    import os
    p.add_argument('--rundir','-C',type=str,default=os.getcwd())
    args=parser.parse_args()

    
    pype1.run(rundir=args.rundir) 
    
    for x in '1PGB 1AKI'.split():
        PDB_ID.put(x)
        # print(PDB_ID.built())
        pype2.run(rundir=args.rundir, target_dir=PDB_ID.value + '/' + str(md_steps))

    pype2.pprint_stats()
main()
