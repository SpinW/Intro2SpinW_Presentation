%% Internal data structure of SpinW
% Import the LiCrO2 crystal structure:
% - spgr: R-3m
% - rhombohedral structure
% - hexagonal settings
% What is the density? spinw.formula()
% Where does all the data comes from? edit atom.dat
%

licro = spinw('licro2.cif');

%% check out internal data structure of SpinW
% What do you see?

properties(spinw)

licro.lattice


%% Test cordinate transformations
% Convert (1,1,0) in the reciprocal lattice of LiCrO2 to A^-1
%
% How are Trl and Tl related?
%

Trl = licro.rl;
Q = [0 1 0];
QinvA = Q*Trl

Tl = licro.basisvector;
pos = [1;0;0];
RA = Tl*pos


%% Rotate a vector
% Rotation angle is given in degree here

R = sw_rotmatd([0 0 1],60);

R*pos

%% Magnetic atoms
% Can you find the magnetic atom positions on the plot and connect them to
% the list in mAtom and in licro.unit_cell?
%

plot(licro)
mAtom = licro.matom

%% Space groups
% Check out symmetry.dat:
% What does it mean P0~=P1?

edit symmetry.dat

%% Bond symmetry
% Generate the bond list of LiCrO2
% Are the bonds generated usign the space group?
% What does it mean?


licro.gencoupling

% The first neighbor bonds

licro.table('bond')


%% Lets do something crazy?
% Lets convert the hexagonal cell to the primitive cell and see what
% happens
%
% Why copy()?
% Check spinw.formula() to see how much the new cell is smaller
%


licroR = copy(licro);
T = licroR.newcell({[-1 1 1]/3 [2 1 1]/3 [-1 -2 1]/3})
plot(licroR,'cellMode','outside')

%% Lets add some matrix
% What do the different inputs do?
% Check licro.matrix
%

licro.addmatrix('label','J','value',1)
licro.addmatrix('label','D','value',[0 1 0])
licro.addmatrix('label','A','value',diag([1 0.5 0.5]))

%% Lets add some matrices and plot
% Do both for the primitive and the hexagonal cell
%
% Are they the same?
%

licroR.gencoupling
licroR.addmatrix('label','J','value',1)
licroR.addcoupling('bond',1,'mat','J')
plot(licroR,'range',[3 3 3],'unit','lu','atomMode','mag','cellMode','inside','atomLegend',false)

%% Plot a single triangular layer from LiCrO2

licro.addcoupling('bond',1,'mat','J')
plot(licro,'range',[15 15 5],'unit','xyz','atomMode','mag','cellMode','none','atomLegend',false)


%% Single ion properties
% 

licro.addaniso('A')
plot(licro)


%% Magnetic structure
% Lets start with a simpler lattice: square lattice
%
% Define a (1/2,1/2,0) magnetic structure using complex magnetization
% vectors
%
% What does [1-1i;1+1i;0] vector means?
%
% How is this stored in spinw object? check spinw.mag_str
%
% Let's create a magnetic supercell! 2x2x1 cell
% what will be the k-vector now? (1/2,1/2,0) or (1,1,0)?
%
% To check that we are doing the right thing, we will keep track of the
% enrgy of the system, after adding a J=1 meV Heisenberg exchange on the
% square lattice.
% Quick Hamiltonian: square.quickham(1)
%
% What was the energy/spin for the original magnetic structure? S=1, J = 1
% meV
%
%


square = sw_model('squareAF');
plot(square,'range',[2 2 1/2])

square.genmagstr('mode','direct','S',[1-1i;1+1i;0],'k',[1/2 1/2 0],'nExt',[1 1 1])
plot(square,'range',[2 2 1/2])

square.quickham(1)

square.energy

% create supercell, which k-vector is right? k is always in the units of
% the crystal reciprocal lattice!

square.genmagstr('mode','direct','S',[1 -1 -1 1;0 0 0 0;0 0 0 0],'nExt',[2 2 1],'k',[0 0 0])
square.energy

square.genmagstr('mode','direct','S',[1 -1 -1 1;0 0 0 0;0 0 0 0],'nExt',[2 2 1],'k',[1/2 1/2 0])
square.energy

%% rotating frame representation
% Let's generate a spiral using the J1-J2 model of a zig-zag chain.
% Check the spinw.mag_str and spinw.magstr and compare them.

ch = spinw;
ch.genlattice('lat_const',[3 5 4])
ch.addatom('r',[0 .3 0],'S',1)
ch.addatom('r',[1/2 0.7 0],'S',1)
ch.quickham([1 1/2])

ch.optmagk
ch.optmagsteep
ch.energy

%%
plot(ch,'range',[5 1 1/2])

%% animation
for ii = 1:100
    ch.genmagstr('mode','rotate','n',[0 0 1],'phid',1)
    swplot.plotmag
    drawnow
end

%%  check spin wave dispersion

spec = ch.spinwave({[0 0 0] [1 0 0] 501});
spec = sw_egrid(spec);
spec = sw_instrument(spec,'dE',0.1,'norm',true);
figure
sw_plotspec(spec,'mode','color')
sw_plotspec(spec,'mode','disp','imag',true)

%%
figure
spec = sw_omegasum(spec,'zeroint',1e-8);
sw_plotspec(spec,'mode','int','axLim',[0 200])

%% rotating frame representation fails
% two counterrotating spirals

ch = spinw;
ch.genlattice('lat_const',[3 7 4])
ch.addatom('r',[1/2 1/4 0],'S',1)
ch.addatom('r',[1/2 3/4 0],'S',1)
ch.addmatrix('label','DM1','value',[0  1 0])
ch.addmatrix('label','DM2','value',[0 -1 0])

ch.gencoupling
ch.addcoupling('mat','DM1','bond',1,'subIdx',1)
ch.addcoupling('mat','DM2','bond',1,'subIdx',2)

optRes = ch.optmagk;
plot(ch,'range',[5 1 1/2])

% rotating frame representation fails!
ch.magstr

% lets check how the structure would look like in a supercell
ch.magstr('nExt',[4 1 1])

% let's generate this supercell and keep
ch.genmagstr('mode','helical','nExt',[4 1 1])

% energy
ch.energy

