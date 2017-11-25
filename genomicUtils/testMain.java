package genomicUtils;

import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;

public class testMain {

	public static void main(String[] args) {

		RegionVector testVector = new RegionVector();

		int realData = 0;

		if (realData == 1) {
			String testString = "TTCGGGCCCAAACACGTGTTATATTATCGAGGTGGAAAACAGCCATATGTATTAGGAGGGGAAACTGAGGCAGAC";
			Region r1 = new Region(50016512, 50016730);
			Region r2 = new Region(50028714, 50028830);
			Region r3 = new Region(50029267, 50029590);

			testVector.addRegion(r1);
			testVector.addRegion(r2);
			testVector.addRegion(r3);
			RegionVector genregions = getTheShitDone(testString, testVector, 585);

			System.out.println(genregions);

		} else {
			String testString = "12345";

			Region r4 = new Region(10, 20); // 0:9
			Region r5 = new Region(30, 40); // 10:19
			Region r6 = new Region(50, 60); // 20:35
			Region r7 = new Region(70, 80); // 0:9
			Region r8 = new Region(100, 140); // 10:19
			Region r9 = new Region(150, 160); // 20:35
			Region r10 = new Region(210, 220); // 0:9
			Region r11 = new Region(230, 240); // 10:19
			Region r12 = new Region(250, 260); // 20:35
			Region r13 = new Region(310, 320); // 0:9
			Region r14 = new Region(330, 340); // 10:19
			Region r15 = new Region(350, 360); // 20:35
			Region r16 = new Region(410, 420); // 0:9
			Region r17 = new Region(430, 440); // 10:19
			Region r18 = new Region(450, 460); // 20:35
			Region r19 = new Region(510, 520); // 0:9
			Region r20 = new Region(530, 540); // 10:19
			Region r21 = new Region(550, 560); // 20:35

			testVector.addRegion(r4);
			testVector.addRegion(r5);
			testVector.addRegion(r6);
			testVector.addRegion(r7);
			testVector.addRegion(r8);
			testVector.addRegion(r9);
			testVector.addRegion(r10);
			testVector.addRegion(r11);
			testVector.addRegion(r12);
			testVector.addRegion(r13);
			testVector.addRegion(r14);
			testVector.addRegion(r15);
			testVector.addRegion(r16);
			testVector.addRegion(r17);
			testVector.addRegion(r18);
			testVector.addRegion(r19);
			testVector.addRegion(r20);
			testVector.addRegion(r21);

			System.out.println(testVector.getRegionLength() + testVector.getSize());

			// RegionVector genregions = getTheShitDone(testString, testVector, 1);
			// System.out.println(genregions);

			for (int i = 1; i < 225; i++) {
				RegionVector genregions = getTheShitDone(testString, testVector, i);
				System.out.println(i + " " + genregions);
			}

		}
	}

	
	
	public static void getTranscriptRegion (String trans, String rev, String fw) {
		
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
				while (rdmStart + FL >= distanceTravelled) // save regions inside of fragment
				{
					distanceTravelled += parent.regions.get(i).getLength();
					if (rdmStart + FL > distanceTravelled) // if FALSE: arrived at last region
					{
						genomicRegions.addRegion(parent.regions.get(i));
					}
					i++;
				}				
			} catch (Exception e) {
				throw new RuntimeException("RandomStart:  " + rdmStart + " FragmentL: " + FL + " DT: " + distanceTravelled + "\n" + parent + "\n" + fragment);
			}

			int pos12;
			Region lastRegion = new Region();
			pos12 = parent.regions.get(i - 1).getX2() - (distanceTravelled - (rdmStart + FL));
			lastRegion = new Region(parent.regions.get(i - 1).getX1(), pos12 - 1);

			genomicRegions.addRegion(lastRegion);

		}
		return genomicRegions;
	}
}
