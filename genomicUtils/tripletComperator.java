package genomicUtils;

public class tripletComperator implements Comperator<Triplet<String, String, Integer>> {

	@Override
	public int compare(Triplet<String, String, Integer> o1, Triplet<String, String, Integer> o2) {

		int diff = 1;
		if(o1.getSecond().equals(o2.getSecond())) {
			diff = 0;
		}

		return diff;
	}

	
	
}
