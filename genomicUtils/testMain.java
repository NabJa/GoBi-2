package genomicUtils;

public class testMain {

	public static void main(String[] args) {

		RegionVector testVector = new RegionVector();
		String testString = "12345";

		Region r1 = new Region(10, 20);
		Region r2 = new Region(25, 30);
		Region r3 = new Region(35, 40);
		Region r4 = new Region(45, 50);
		Region r5 = new Region(55, 80);

		testVector.addRegion(r1);
		testVector.addRegion(r2);
		testVector.addRegion(r3);
		testVector.addRegion(r4);
		testVector.addRegion(r5);

		RegionVector genregions = getTheShitDone(testString, testVector, 2);

		System.out.println(genregions);
	}

	public static RegionVector getTheShitDone(String fragment, RegionVector parent, int rdmStart) {
		RegionVector genomicRegions = new RegionVector();

		int FL = fragment.length();
		int distanceTravelled = 0;
		int i = 0;

		while (rdmStart > distanceTravelled) // skip regions if start isnt in it
		{
			distanceTravelled += parent.regions.get(i).getLength();
			i++;
		}

		int pos1;
		if (i - 1 == 0) // if start is in first region
		{
			pos1 = parent.regions.get(i - 1).getX1() + rdmStart;
		} else {
			pos1 = parent.regions.get(i - 1).getX1()
					+ (rdmStart - (distanceTravelled - parent.regions.get(i - 1).getLength()) - 1);
		}

		if (rdmStart + FL < distanceTravelled) // if end ist in first region
		{
			int pos12 = pos1 + FL;
			Region onlyRegion = new Region(pos1, pos12);
			genomicRegions.addRegion(onlyRegion);

		} else {
			Region firstRegion = new Region(pos1, parent.regions.get(i - 1).getX2());
			genomicRegions.addRegion(firstRegion);

			while (rdmStart + FL > distanceTravelled) // save regions inside of fragment
			{
				distanceTravelled += parent.regions.get(i).getLength();
				if (rdmStart + FL > distanceTravelled) // if FALSE: arrived at last region
				{
					genomicRegions.addRegion(parent.regions.get(i));
				}
				i++;
			}
			int pos12 = parent.regions.get(i - 1).getX2() - (distanceTravelled - (rdmStart + FL));
			Region lastRegion = new Region(parent.regions.get(i - 1).getX1(), pos12);
			genomicRegions.addRegion(lastRegion);
		}
		return genomicRegions;
	}
}
