package genomicUtils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;

public class RegionVector {

	public String id; // transcript ID
	public int x1 = Integer.MAX_VALUE; // transcript start
	public int x2 = 0; // transcript end

	public ArrayList<Region> regions = new ArrayList<Region>(); // maybe Vector instead of ArrayList???

	public RegionVector() {
	}

	public RegionVector(String id) {
		this.id = id;
	}

	public RegionVector(String id, int x1, int x2) {
		this.id = id;
		this.x1 = x1;
		this.x2 = x2;
	}

	public String getID() {
		return id;
	}

	public int getX1() {
		return x1;
	}

	public int getX2() {
		return x2;
	}

	public int getLength() {
		return x2 - x1;
	}

	public int getRegionLength() {
		int len = 0;
		for (Region r : regions) {
			len += r.getLength();
		}
		return len;
	}

	/**
	 * adds the Region in a sorted way. Smallest start first.
	 * 
	 * @param region
	 */
	public void addRegion(Region region) {
		if (regions.size() == 0) {
			regions.add(region);
		} else {
			for (int i = 0; i < regions.size(); i++) {
				if (regions.get(i).getX1() < region.getX1())
					continue;
				regions.add(i, region);
				return;
			}
			regions.add(region);
		}

		if (region.getX1() < x1) {
			x1 = region.getX1();
		}
		if (region.getX2() > x2) {
			x2 = region.getX2();
		}
	}

	public int[] getTransLoc() {
		int[] region = { x1, x2 };
		return region;
	}

	/**
	 * Returns an ArrayList<Region> that are the inverse of given RegionVector.
	 * Example: Input{1,2; 3,4; 5,6} -> Output{2,3; 4,5}
	 * 
	 * @return
	 */
	public ArrayList<Region> arrayInverse() {

		ArrayList<Region> introns = new ArrayList<Region>();

		for (int i = 0; i < regions.size() - 1; i++) {
			Region intron = new Region(regions.get(i).getX2(), regions.get(i + 1).getX1());
			introns.add(intron);
		}

		return introns;
	}

	/**
	 * Returns an RegionVector that are the inverse of given RegionVector. Has to be
	 * save in new Variable!!! Example: Input{1,2; 3,4; 5,6} -> Output{2,3; 4,5}
	 * 
	 * @return
	 */
	public RegionVector inverse() {

		RegionVector introns = new RegionVector();

		for (int i = 0; i < regions.size() - 1; i++) {
			Region intron = new Region(regions.get(i).getX2() + 1, regions.get(i + 1).getX1());
			introns.addRegion(intron);
		}

		introns.id = id;

		return introns;
	}

	public void printRegions() {
		System.out.println("ID: " + id);
		for (Region r : regions) {
			System.out.println(r.getX1() + " " + r.getX2() + " ");
		}
		System.out.println();
	}

	public int getSize() {
		int i = 0;
		for (Region r : regions) {
			i++;
		}
		return i;
	}

	@Override
	public int hashCode() {
		return ((x1 * 104723) % 104729) + ((x2 * 104717) % 104711);
	}

	@Override
	public String toString() {
		String result = id;
		for (int i = 0; i < regions.size(); i++) {
			result += " " + regions.get(i);
		}
		return result;
	}
}