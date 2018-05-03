import static java.lang.Math.pow;
import static java.lang.StrictMath.sin;
import static java.lang.StrictMath.sqrt;
import static java.lang.Math.cos;

public class MHD {

    private Vector waveVector;
    private Vector velocity;
    private Vector magneticField;


    public MHD(){}

    public MHD(Vector wave_vector, Vector velocity, Vector magneticField) {
        this.waveVector = wave_vector;
        this.velocity = velocity;
        this.magneticField = magneticField;
    }

    public Vector getWaveVector() {
        return waveVector;
    }

    public void setWaveVector(Vector waveVector) {
        this.waveVector = waveVector;
    }

    public Vector getVelocity() {
        return velocity;
    }

    public void setVelocity(Vector velocity) {
        this.velocity = velocity;
    }

    public Vector getMagneticField() {
        return magneticField;
    }

    public void setMagneticField(Vector magneticField) {
        this.magneticField = magneticField;
    }

    public MHD linearTerms(){
        final double nu_k=0.001;
        final double nu_m=0.001;
        final double omega=-2.0/3.0;
        final double alpha=0.0;

        double Q2 = pow(waveVector.getX(),2)+ pow(waveVector.getY(),2) + pow(waveVector.getZ(),2);
        double Q = sqrt(Q2);
        double Qa = waveVector.getY()*sin(alpha) + waveVector.getZ()*cos(alpha);
        double nx = waveVector.getX()/Q;
        double ny = waveVector.getY()/Q;
        double nz = waveVector.getZ()/Q;


        double cm = 2.0*(ny*velocity.getX() - omega*(nx*velocity.getY()-ny*velocity.getX()));
        double vx =                   cm*nx + 2.0*omega*velocity.getY() + Qa*magneticField.getX() - nu_k*Q2*velocity.getX();
        double vy =-velocity.getX() + cm*ny - 2.0*omega*velocity.getX() + Qa*magneticField.getY() - nu_k*Q2*velocity.getY();
        double vz =                   cm*nz                             + Qa*magneticField.getZ() - nu_k*Q2*velocity.getZ();
        Vector newVelocity = new Vector(vx,vy,vz);

        double bx =                      - Qa*velocity.getX() - nu_m*Q2*magneticField.getX();
        double by = magneticField.getX() - Qa*velocity.getY() - nu_m*Q2*magneticField.getY();
        double bz=                       - Qa*velocity.getZ() - nu_m*Q2*magneticField.getZ();

        Vector newMagneticField = new Vector(bx,by,bz);
        MHD next_mhd = new MHD(waveVector,newVelocity,newMagneticField);
        return next_mhd;
    }
}
