package genomicUtils;

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
			String testString = "1234";

			Region r4 = new Region(10, 20); // 0:9
			Region r5 = new Region(30, 40); // 10:19
			Region r6 = new Region(50, 60); // 20:35

			testVector.addRegion(r4);
			testVector.addRegion(r5);
			testVector.addRegion(r6);

			RegionVector genregions = getTheShitDone(testString, testVector, 26);

			System.out.println(genregions);
		}
	}

	public static RegionVector getTheShitDone(String fragment, RegionVector parent, int rdmStart) {

		RegionVector genomicRegions = new RegionVector();

		int FL = fragment.length();
		int distanceTravelled = 0;
		int i = 0;

		while (rdmStart >= distanceTravelled) // skip regions if start isnt in it
		{
			distanceTravelled += parent.regions.get(i).getLength() + 1;
			i++;
		}

		int pos1 = 0;

		if (i - 1 == 0) // if start is in first region
		{
			pos1 = parent.regions.get(i - 1).getX1() + rdmStart;
		} else if (i == 0) {
			pos1 = parent.regions.get(i).getX1();
		} else {
			pos1 = parent.regions.get(i - 1).getX1() + (rdmStart - (distanceTravelled - parent.regions.get(i - 1).getLength() - i));
		}

		if (rdmStart + FL < distanceTravelled - i ) // if end is in first region
		{
			int pos12 = pos1 + FL;
			Region onlyRegion = new Region(pos1, pos12);
			genomicRegions.addRegion(onlyRegion);

		} else {
			Region firstRegion = new Region(pos1, parent.regions.get(i - 1).getX2());
			genomicRegions.addRegion(firstRegion);

			while (rdmStart + FL > distanceTravelled + i) // save regions inside of fragment
			{
				distanceTravelled += parent.regions.get(i).getLength();
				if (rdmStart + FL > distanceTravelled) // if FALSE: arrived at last region
				{
					genomicRegions.addRegion(parent.regions.get(i));
				}
				i++;
			}
			if (i != parent.regions.size()) {

				int pos12;
				Region lastRegion = new Region();
				if (distanceTravelled - i != (rdmStart + FL)) {
					pos12 = parent.regions.get(i - 1).getX2() - (distanceTravelled - (rdmStart + FL));
					lastRegion = new Region(parent.regions.get(i - 1).getX1(), pos12 + 1);
				} else {
					pos12 = parent.regions.get(i).getX2() - distanceTravelled + i;
					lastRegion = new Region(parent.regions.get(i).getX1(), pos12 + 1);
				}
				genomicRegions.addRegion(lastRegion);
			}
		}
		return genomicRegions;
	}
}
