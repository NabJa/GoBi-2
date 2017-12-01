package genomicUtils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

public class MutationDistribution {
	String outputDestination;
	FileWriter file;
	BufferedWriter writer;

	public MutationDistribution(File outputDestination) {

		try {
			this.file = new FileWriter(outputDestination);
			this.writer = new BufferedWriter(file);
		} catch (Exception e) {
			throw new RuntimeException("got error while creating mutationDis", e);
		}

		try {
			writer.write(
					"Mutation_Pos");
			writer.newLine();
		} catch (Exception e) {
			throw new RuntimeException("Got error while printing mutheader header", e);
		}

	}

	public void writeMutOut(ArrayList<Integer> fw_mut, ArrayList<Integer> rw_mut) {

		for (int i = 0; i < fw_mut.size(); i++) {
			try {
				writer.write(fw_mut.get(i));
				writer.newLine();
			} catch (Exception e) {
				throw new RuntimeException("Error while printing mutation distribution at: " + i + "\n" + fw_mut, e);
			}
		}

		for (int i = 0; i < rw_mut.size(); i++) {

			try {
				writer.write(rw_mut.get(i));
				writer.newLine();
			} catch (Exception e) {
				throw new RuntimeException("Error while printing mutation distribution at: " + i + "\n" + rw_mut, e);
			}
		}
	}

	public void close() {
		try {
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

}
