package algorithm;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.ArrayList;

public class DPAlgorithm {
	int maxTADbits = 10000000;
	int maxTADLength;
	int windowBits = 10000000;
	int windowSize;
	RealMatrix graph;
	RealMatrix rowSum;
	int numOfVertices;
	int dim;
	double sumOfDegrees;
	ArrayList<Integer> zeroRowList;
	int binSize;
	
	public DPAlgorithm(RealMatrix mtx, int binSize, int windowBits, int maxTADbits) {
		this.binSize = binSize;
		this.windowBits = windowBits;
		this.windowSize = windowBits / binSize;
		this.maxTADbits = maxTADbits;
		maxTADLength = maxTADbits / binSize;
		graph = mtx.copy();
		mtx = null;
		dim = graph.getRowDimension();
		for(int i=0; i<dim; i++)
			graph.setEntry(i, i, 0);
		numOfVertices = dim;
		RealMatrix ones = MatrixUtils.createRealMatrix(dim, 1);
		for(int i=0; i<dim; i++)
			ones.setEntry(i, 0, 1);
		rowSum = graph.multiply(ones);
		sumOfDegrees = ones.transpose().multiply(rowSum).getEntry(0, 0);
		zeroRowList = new ArrayList<>();
		for(int i=0; i<dim; i++)
			if(rowSum.getEntry(i, 0) == 0)
				zeroRowList.add(i);
		for(int i=0; i<dim; i++) {
			int nonZeroRows = dim - zeroRowList.size();
			if(rowSum.getEntry(i, 0) == 0) {
				if(sumOfDegrees > 0) {
					graph.setEntry(i, i, sumOfDegrees/nonZeroRows);		//set the diagonal of zero rows a number to make each row have an average number.
					rowSum.setEntry(i, 0, sumOfDegrees/nonZeroRows);
				}
			}
		}
	}
	
	/*
	 * get structure entropy of a node on the tree.
	 * */
	public double getNodeEntro(double cut, double volume, double parVolume) {
		double result = - (cut/sumOfDegrees) * (Math.log(volume/parVolume)/Math.log(2));
		return result;
	}
	
	/*
	 * deDoc2.w: minimize 2D structural entropy using dynamic programming algorithm and get TLD partition.
	 * */
	public ArrayList<ArrayList<Integer>> getPartition2D(){	
		RealMatrix volumeInsideMtx = new BlockRealMatrix(dim, dim);
		for(int i=0; i<dim; i++)
			volumeInsideMtx.setEntry(i, i, graph.getEntry(i, i));
		for(int i=0; i<dim; i++) {
			for(int k=1; i+k<dim && k<=maxTADLength; k++) {
				int j = i+k;
				double valueij = volumeInsideMtx.getEntry(i, j-1);
				for(int m=i; m<j; m++) {
					valueij += graph.getEntry(j, m);
					valueij += graph.getEntry(m, j);
				}
				valueij += graph.getEntry(j, j);
				volumeInsideMtx.setEntry(i, j, valueij);
			}
		}
		
		RealMatrix cutMtx = new BlockRealMatrix(dim, dim);
		for(int i=0; i<dim; i++) {
			for(int j=i; j<dim && j-i<=maxTADLength; j++) {
				double cutij = 0;
				for(int k=i; k<=j; k++)
					cutij += rowSum.getEntry(k, 0);
				cutij -= volumeInsideMtx.getEntry(i, j);
				cutMtx.setEntry(i, j, cutij);
			}
		}
		
		RealMatrix communityEntro = new BlockRealMatrix(dim, dim);
		for(int i=0; i<dim; i++) {
			for(int j=i; j<dim && j-i<=maxTADLength; j++) {
				double entroij = 0;
				double volumeij = 0;
				for(int k=i; k<=j; k++)
					volumeij += rowSum.getEntry(k, 0);
				if(volumeij > 0) {
					entroij += getNodeEntro(cutMtx.getEntry(i, j), volumeij, sumOfDegrees);
					for(int k=i; k<=j; k++) {
						entroij += getNodeEntro(cutMtx.getEntry(k, k), rowSum.getEntry(k, 0), volumeij);
					}
				}
				else {
					entroij = 0;
				}
				communityEntro.setEntry(i, j, entroij);
			}
		}
		
		ArrayList<Double> strucEntroList = new ArrayList<>();
		for(int i=0; i<=dim+1; i++)
			strucEntroList.add(-1.0);
		ArrayList<Integer> decisionList = new ArrayList<>();
		for(int i=0; i<=dim+1; i++)
			decisionList.add(0);
		strucEntroList.set(dim, 0.0);
		for(int i=dim-1; i>=0; i--) {
			double min=Double.MAX_VALUE;
			int minIndex = -1;
			for(int k=0; k<dim-i && k<=maxTADLength; k++) {
				double strucEntroik = communityEntro.getEntry(i, i+k) + strucEntroList.get(i+k+1);
				if(strucEntroik < min) {
					min = strucEntroik;
					minIndex = k;
				}
			}
			if(minIndex == -1)
				throw new RuntimeException("no right TAD boundary found for bin " + i);
			strucEntroList.set(i, min);
			decisionList.set(i, i + minIndex + 1);
		}
		
		ArrayList<Integer> boundaries = new ArrayList<>();
		int boundaryi = 0;
		while(boundaryi<dim) {
			boundaries.add(boundaryi);
			boundaryi = decisionList.get(boundaryi);
		}
		ArrayList<ArrayList<Integer>> partition2D = boundariesToPartition(boundaries);
		
		return partition2D;
	}
	
	
	/*
	 * deDoc2.s: split Hi-C contact matrix into proper windows and minimize 2D structural entropy of the graph inside the windows using dynamic algorithm.
	 * this function predict TLDs in one window each time.
	 * */
	public ArrayList<ArrayList<Integer>> getPartitionWindow(int nextWindowStart){
		if(dim - nextWindowStart < 2*windowSize) {
			RealMatrix windowMtx = graph.getSubMatrix(nextWindowStart, dim-1, nextWindowStart, dim-1);
			ArrayList<ArrayList<Integer>> partitionWindow = new DPAlgorithm(windowMtx, binSize, windowBits, maxTADbits).getPartition2D();
			for(ArrayList<Integer> community : partitionWindow) {
				for(int j=0; j<community.size(); j++) {
					community.set(j, community.get(j)+nextWindowStart);
				}
			}
			return partitionWindow;
		}
		else {
			RealMatrix windowMtx = graph.getSubMatrix(nextWindowStart, nextWindowStart+windowSize-1, nextWindowStart, nextWindowStart+windowSize-1);
			ArrayList<ArrayList<Integer>> partitionWindow = new DPAlgorithm(windowMtx, binSize, windowBits, maxTADbits).getPartition2D();
			for(ArrayList<Integer> community : partitionWindow) {
				for(int j=0; j<community.size(); j++) {
					community.set(j, community.get(j)+nextWindowStart);
				}
			}
			nextWindowStart = partitionWindow.get(partitionWindow.size()-1).get(0);
			ArrayList<ArrayList<Integer>> partitionLeft = getPartitionWindow(nextWindowStart);
			ArrayList<ArrayList<Integer>> result = new ArrayList<>();
			for(int i=0; i<partitionWindow.size()-1; i++)
				result.add(partitionWindow.get(i));
			for(int i=0; i<partitionLeft.size(); i++)
				result.add(partitionLeft.get(i));
			return result;
		}
	}
	
	/*
	 * transform boundary format to partition format.
	 * */
	public ArrayList<ArrayList<Integer>> boundariesToPartition(ArrayList<Integer> boundaries){
		ArrayList<ArrayList<Integer>> partition = new ArrayList<>();
		
		int pos = 0;
		for(int i=1; i<boundaries.size(); i++) {
			int boundaryi = boundaries.get(i);
			ArrayList<Integer> row = new ArrayList<>();
			for(; pos<boundaryi; pos++) {
				row.add(pos);
			}
			partition.add(row);
		}
		ArrayList<Integer> row = new ArrayList<>();
		for(; pos<dim; pos++) {
			row.add(pos);
		}
		partition.add(row);
		partition = trimZeros(partition);
		return partition;
	}
	
	
	/*
	 * deal with the gap regions.
	 * */
	public ArrayList<ArrayList<Integer>> trimZeros(ArrayList<ArrayList<Integer>> partitionOrigin){
		ArrayList<ArrayList<Integer>> partition = new ArrayList<>();
		for(int i=0; i<partitionOrigin.size(); i++) {
			ArrayList<Integer> communityi = partitionOrigin.get(i);
			if(communityi.size()>=3) {
				boolean breakCommunity = true;
				for(int j=1; j<communityi.size()-1; j++) {
					if(!zeroRowList.contains(communityi.get(j))) {
						breakCommunity = false;
					}
				}
				if(breakCommunity) {
//					System.out.println("break community from "+communityi.get(0)+" to "+communityi.get(communityi.size()-1));
					for(int j=0; j<communityi.size(); j++) {
						ArrayList<Integer> com = new ArrayList<>();
						com.add(communityi.get(j));
						partition.add(com);
					}
					continue;
				}
			}
			int nonZeroStart = -1;
			for(int j=0; j<communityi.size(); j++) {
				if(!zeroRowList.contains(communityi.get(j))) {
					nonZeroStart = j;
					break;
				}
			}
			int nonZeroEnd = -1;
			for(int j=communityi.size()-1; j>=0; j--) {
				if(!zeroRowList.contains(communityi.get(j))) {
					nonZeroEnd = j;
					break;
				}
			}
			if(nonZeroStart == -1 || nonZeroEnd == -1) {
				for(int j=0; j<communityi.size(); j++) {
					ArrayList<Integer> com = new ArrayList<>();
					com.add(communityi.get(j));
					partition.add(com);
				}
			}
			else {
				for(int j=0; j<nonZeroStart; j++) {
					ArrayList<Integer> com = new ArrayList<>();
					com.add(communityi.get(j));
					partition.add(com);
				}
				ArrayList<Integer> nonZeroCom = new ArrayList<>();
				for(int j=nonZeroStart; j<=nonZeroEnd; j++) {
					nonZeroCom.add(communityi.get(j));
				}
				partition.add(nonZeroCom);
				for(int j=nonZeroEnd+1; j<communityi.size(); j++) {
					ArrayList<Integer> com = new ArrayList<>();
					com.add(communityi.get(j));
					partition.add(com);
				}
			}
		}
		return partition;
	}
	
	/*
	 * calculate 2D structural entropy of a graph given a partition.
	 * */
	public double calStrucEntro2D(ArrayList<ArrayList<Integer>> partition) {
		double L = 0;
		double D = 0;
		double V = 0;
		double G = 0;
		for(int i=0;i<partition.size();i++) {
			ArrayList<Integer> communityi = partition.get(i);
			double volumnInside = 0;
			double gi = 0;
			double vi = 0;
			double selfloopi = 0;
			int comStart = communityi.get(0);
			int comEnd = communityi.get(communityi.size()-1);
			for(int j=comStart;j<=comEnd;j++) {
				for(int k=comStart;k<=comEnd;k++) {
					volumnInside += graph.getEntry(j, k);
				}
				gi += rowSum.getEntry(j, 0);
				vi += rowSum.getEntry(j, 0);
				selfloopi += graph.getEntry(j, j);
				double dj = rowSum.getEntry(j, 0);
				if(dj>0)
					D += (dj-graph.getEntry(j, j))*(Math.log(dj)/Math.log(2));
			}
			gi -= volumnInside;
			if(vi>0) {
				V += (vi-selfloopi)*(Math.log(vi)/Math.log(2));
				G += gi*(Math.log(vi/sumOfDegrees)/Math.log(2));
			}
		}
		L = (V-D-G)/sumOfDegrees;
		return L;
	}
	
	
	/*
	 * calculate the 1D structural entropy of the graph.
	 * */
	public double calStrucEntro1D() {
		double entropy = 0;
		for(int i=0;i<numOfVertices;i++) {
			double di = rowSum.getEntry(i, 0);
			if(di>0)
				entropy -= (di-graph.getEntry(i, i))*(Math.log(di/sumOfDegrees)/Math.log(2));
		}
		entropy /= sumOfDegrees;
		return entropy;
	}
	
	/*
	 * write the TLDs result to files.
	 * */
	public void putPartition(ArrayList<ArrayList<Integer>> partition, String fileName) {
		try {
			BufferedWriter bw = new BufferedWriter(new FileWriter(fileName));
			for(int i=0;i<partition.size();i++) {
				ArrayList<Integer> communityi = partition.get(i);
				for(int j=0;j<communityi.size();j++) {
					bw.write((communityi.get(j)+1) + " ");
				}
				bw.write("\n");
			}
			bw.close();
		}catch(Exception e) {
			e.printStackTrace();
		}
	}
	
	
	/*
	 * calculate deDoc2.w and deDoc2.s and write the TLDs results.
	 * */
	public void putTLD(String fileName) {
		ArrayList<ArrayList<Integer>> partition2D = getPartition2D();
		putPartition(partition2D, fileName + ".TAD");
		ArrayList<ArrayList<Integer>> partitionWindow = getPartitionWindow(0);
		putPartition(partitionWindow, fileName+".window.TAD");
	}
}
