## qCom.sh is a command which sends each commandline to our local grid

## directory to store results files

## what happens when we vary parameters one at a time?
MAF='0.03,0.05,0.1,0.15,0.2,0.3,0.4'
P1='5e-6,1e-6,5e-7,1e-7,5e-8,1e-8,5e-9,1e-9,1e-10,1e-12,1e-15,1e-20'
P2='1e-3,5e-4,2e-4,1e-4,5e-5,2e-5,1e-5,5e-6,2e-6,1e-6,5e-7,2e-7,1e-7,5e-8,2e-8,1e-8'
WTCCC='N0=3000 N1=2000 N2=2000'
ICHIP='N0=9000 N1=6000 N2=7000'
N='2e3,3e3,4e3,5e3,6e3,7e3,8e3,9e3,10e3'

qCom.sh ./simulations.R --args $WTCCC q0=$MAF p1=5e-8 p2=1e-4 file=wtccc-maf
qCom.sh ./simulations.R --args $ICHIP q0=$MAF p1=5e-8 p2=1e-4 file=ichip-maf
qCom.sh ./simulations.R --args $WTCCC q0=0.2 p1=$P1 p2=1e-4 file=wtccc-p1
qCom.sh ./simulations.R --args $ICHIP q0=0.2 p1=$P1 p2=1e-4 file=ichip-p1
qCom.sh ./simulations.R --args $WTCCC q0=0.2 p1=5e-8 p2=$P2 file=wtccc-p2
qCom.sh ./simulations.R --args $ICHIP q0=0.2 p1=5e-8 p2=$P2 file=ichip-p2

qCom.sh ./simulations.R --args N0=$N $WTCCC q0=0.2 p1=5e-8 p2=1e-4 file=wtccc-N0
qCom.sh ./simulations.R --args N1=$N $WTCCC q0=0.2 p1=5e-8 p2=1e-4 file=wtccc-N1
qCom.sh ./simulations.R --args N2=$N $WTCCC q0=0.2 p1=5e-8 p2=1e-4 file=wtccc-N2
qCom.sh ./simulations.R --args N0=$N $ICHIP q0=0.2 p1=5e-8 p2=1e-4 file=ichip-N0
qCom.sh ./simulations.R --args N1=$N $ICHIP q0=0.2 p1=5e-8 p2=1e-4 file=ichip-N1
qCom.sh ./simulations.R --args N2=$N $ICHIP q0=0.2 p1=5e-8 p2=1e-4 file=ichip-N2
