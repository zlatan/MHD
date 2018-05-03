import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static java.util.stream.Collectors.toList;

public class Disk {
    public static void main(String[] args) {


        List<MHD> mhds = new ArrayList<MHD>();

        List<Double> qx = IntStream
                .rangeClosed(-5, 5)
                .mapToDouble(i -> i)
                .map(x -> x / 10.0)
                .filter(x -> x != 0.0)
                .boxed()
                .collect(toList());

        List<Double> qy = IntStream
                .rangeClosed(-5, 5)
                .mapToDouble(i -> i)
                .map(x -> x / 10.0)
                .filter(x -> x != 0.0)
                .boxed()
                .collect(toList());

        List<Double> qz = IntStream
                .rangeClosed(-5, 5)
                .mapToDouble(i -> i)
                .map(x -> x / 10.0)
                .filter(x -> x != 0.0)
                .boxed()
                .collect(toList());


        qx.forEach( x -> {
            qy.forEach( y -> {
                qz.forEach( z -> {
                    Vector waveVector = new Vector(x, y, z);
                    Vector velocity = new Vector(1.0, 1.0, 0.0);
                    Vector magneticField = new Vector(1.0, 1.0, 0.0);
                    mhds.add(new MHD(waveVector, velocity, magneticField));
                });
            });
        });

        List<MHD> firstPass = mhds.stream()
                                .map(MHD::linearTerms)
                                .collect(Collectors.toList());

   }
}
