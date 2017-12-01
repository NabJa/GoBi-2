package genomicUtils;

import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;

import genomeExtractor.GenomeSequenceExtractor;

public class testMain {

	public static void main(String[] args) {

		RegionVector testVector = new RegionVector();

		
		Region r3 = new Region(764,774);
		Region r1 = new Region(905, 998);
		Region r2 = new Region(812, 905);
		Region r5 = new Region(812, 905);
		Region r4 = new Region(774,867);
		
		testVector.addRegion(r3);
		testVector.addRegion(r2);
		testVector.addRegion(r1);		
		testVector.addRegion(r4);
		testVector.addRegion(r5);
		
		
		System.out.println(testVector);
		
		String head = "ENST00000523865 cdna:known chromosome:GRCh37:22:24108061:24108651:-1 gene:ENSG00000250479 gene_biotype:protein_coding transcript_biotype:retained_intron";
		boolean start = head.startsWith(">");
		System.out.println(start);
		
	}

	public static RegionVector getTheShitDone(String fragment, RegionVector parent, int rdmStart) {
		RegionVector genomicRegions = new RegionVector();

		int FL = fragment.length();
		int distanceTravelled = 0;
		int i = 0;

		while (rdmStart > distanceTravelled) // skip regions if start isnt in it
		{
			distanceTravelled += parent.regions.get(i).getLength() + 1;
			i++;
		}

		int pos1 = 0;

		if (i - 1 == 0) // if start is in first region
		{
			pos1 = parent.regions.get(i - 1).getX1() + rdmStart - 1;
		} else if (i == 0) {
			pos1 = parent.regions.get(i).getX1();
		} else {
			pos1 = parent.regions.get(i - 1).getX1()
					+ (rdmStart - (distanceTravelled - parent.regions.get(i - 1).getLength()));
		}

		if (rdmStart + FL <= distanceTravelled + 1) // if end is in first region
		{
			int pos12 = pos1 + FL;
			Region onlyRegion = new Region(pos1, pos12);
			genomicRegions.addRegion(onlyRegion);

		} else {
			Region firstRegion = new Region(pos1, parent.regions.get(i - 1).getX2() + 1);
			genomicRegions.addRegion(firstRegion);

			try {
				while (rdmStart + FL > distanceTravelled) // save regions inside of fragment
				{
					distanceTravelled += parent.regions.get(i).getLength();
					if (rdmStart + FL > distanceTravelled) // if FALSE: arrived at last region
					{
						Region middleRegion = new Region(parent.regions.get(i).getX1(), parent.regions.get(i).getX2() + 1);
						genomicRegions.addRegion(middleRegion);
					}
					i++;
				}
			} catch (Exception e) {
				throw new RuntimeException("RandomStart:  " + rdmStart + " FragmentL: " + FL + " DT: "
						+ distanceTravelled + "\n" + parent + "\n" + fragment);
			}

			if (genomicRegions.getRegionLength() - 1 != FL) {

				int pos12;
				Region lastRegion = new Region();
				pos12 = parent.regions.get(i - 1).getX2() - (distanceTravelled - (rdmStart + FL));
				lastRegion = new Region(parent.regions.get(i - 1).getX1(), pos12 - 1);

				genomicRegions.addRegion(lastRegion);
			}
		}
		return genomicRegions;
	}
}
