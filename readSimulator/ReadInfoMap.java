package readSimulator;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import genomicUtils.Region;
import genomicUtils.RegionVector;

public class ReadInfoMap {

	String outputDestination;
	FileWriter file;
	BufferedWriter writer;

	public ReadInfoMap(String outputDestination) {

		this.outputDestination = outputDestination;
		try {
			this.file = new FileWriter(outputDestination);
			this.writer = new BufferedWriter(file);
		} catch (Exception e) {
			throw new RuntimeException("got error while creating Mapinfo", e);
		}
		
		try {
			writer.write("readid\tchr_id\tgene_id\ttranscript_id\tt_fw_regvec\tt_rw_regvec\tfw_regvec\trw_regvec\tfw_mut\trw_mut");
			writer.newLine();
		} catch (Exception e) {
			throw new RuntimeException("Got error while printing mapinfo header" , e);
		}

	}
	
	public void writeMapinfo(long id, String chr, String gene, String transcript, RegionVector fw_regvec, RegionVector rw_regvec, Region t_fw_regvec,
			Region t_rw_regvec, ArrayList<Integer> fw_mut, ArrayList<Integer> rw_mut)	{
		
		String fwRegvec =  fw_regvec.getX1() + "-" + fw_regvec.getX2();
		String rwRegvec =  rw_regvec.getX1() + "-" + rw_regvec.getX2();
		String t_fwRegvec =  t_fw_regvec.getX1() + "-" + t_fw_regvec.getX2();
		String t_rwRegvec =  t_rw_regvec.getX1() + "-" + t_rw_regvec.getX2();

		String fwMut = "";
		for (int fmut : fw_mut) {
			fwMut += ("," + fmut);
		}
		if(fwMut.length() > 1) {
			fwMut = fwMut.substring(1, fwMut.length());			
		}
		
		String rwMut = "";
		for (int rmut : rw_mut) {
			rwMut += ("," + rmut);
		}
		if(rwMut.length() > 1) {
			rwMut = rwMut.substring(1, rwMut.length());			
		}
		
		try {
			
			writer.write(id + "\t" + chr  + "\t" + gene + "\t" + transcript + "\t" + fwRegvec +"\t" + rwRegvec + "\t" +  t_fwRegvec +"\t" + t_rwRegvec + "\t" + fwMut + "\t" + rwMut );
			writer.newLine();
			
		} catch (Exception e) {
			throw new RuntimeException("Error while printing Readinfo Map: ", e);
		}
		
	}
	
	public void closeMapinfo () {
		try {
			writer.close();
		} catch (Exception e) {
			throw new RuntimeException("Got error while closing mapinfo writer", e);
		}
	}

}
