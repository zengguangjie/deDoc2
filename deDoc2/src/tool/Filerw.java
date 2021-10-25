package tool;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.Scanner;

import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.RealMatrix;


public class Filerw {
	/*
	 * get matrix object from matrix file, entries should be tab split.
	 * */
	public static RealMatrix getMatrix(String fileName) {
		RealMatrix mtx = null;
		try {
			Scanner scan = new Scanner(new BufferedReader(new InputStreamReader(new FileInputStream(fileName), "UTF-8")));
			int i = 0;
			while(scan.hasNextLine()) {
				String[] line = scan.nextLine().trim().split("\t");
				if(mtx == null)
					mtx = new BlockRealMatrix(line.length, line.length);
				for(int j=0;j<line.length;j++) {
					double valueij = Double.parseDouble(line[j]);
					mtx.setEntry(i, j, valueij);
				}
				i++;
			}
		}catch(Exception e) {
			e.printStackTrace();
		}
		if(mtx == null)
			System.err.println(fileName + " is empty!");
		return mtx;
	}
	
	/*
	 * get matrix object from sparse Hi-C format file, the format is the same as deDoc.
	 * */
	public static RealMatrix getMatrix_sparse(String fileName) {
		RealMatrix mtx = null;
		try {
			Scanner scan = new Scanner(new BufferedReader(new InputStreamReader(new FileInputStream(fileName), "UTF-8")));
			if(scan.hasNextLine()) {
				int dim = Integer.parseInt(scan.nextLine().trim());
				mtx = new BlockRealMatrix(dim, dim);
			}
			while(scan.hasNextLine()) {
				String[] line = scan.nextLine().trim().split(" ");
				int i = Integer.parseInt(line[0]) - 1;		//node id start from 1.
				int j = Integer.parseInt(line[1]) - 1;
				double weight = Double.parseDouble(line[2]);
				if(mtx.getEntry(i, j) == 0) {
					mtx.setEntry(i, j, weight);
					mtx.setEntry(j, i, weight);
				}
				else {
					if(Math.abs(mtx.getEntry(i, j)-weight) >= 1e-6) {
						System.err.println("Input Hi-C matrix is not symmetric matrix, the first read will be used. The wrong row is "+i+" "+j+" "+weight);
						mtx.setEntry(i, j, weight);
						mtx.setEntry(j, i, weight);
					}
				}
			}
		}catch(Exception e) {
			e.printStackTrace();
		}
		if(mtx == null)
			System.err.println(fileName + " is empty!");
		return mtx;
	}
	
	/*
	 * write matrix object to matrix file.
	 * */
	public static void putMatrix(RealMatrix mtx, String fileName) {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(fileName));
			for(int i=0;i<mtx.getRowDimension();i++) {
				for(int j=0;j<mtx.getColumnDimension()-1;j++) {
					bw.write(mtx.getEntry(i, j) + "\t");
				}
				bw.write(mtx.getEntry(i, mtx.getColumnDimension()-1)+"\n");
			}
			bw.close();
		}catch(Exception e) {
			e.printStackTrace();
		}
	}
	
}
