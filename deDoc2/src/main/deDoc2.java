package main;


import tool.Filerw;
import tool.Impute;
import algorithm.DPAlgorithm;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.math3.linear.RealMatrix;

public class deDoc2 {
	static String inputFileName;
	static boolean sparseformat = false;
	static int binsize;
	static String outputFileName;
	static int windowsize = 10000000;
	static int maxTLDsize = 10000000;
	static boolean rwr = false;
	static double rp = 0.5;
	static double tol = 0.01;
	static boolean NDI = false;
	public static void main(String[] args) throws Exception {
		CommandLineParser parser = new DefaultParser();
		
		Options options = new Options();
		Option inputfileopt = Option.builder("inputfile").argName("file").hasArg().required(true).desc("input Hi-C data file name").build();
		Option sparseformatopt = new Option("sparseformat", "input file format is sparse format");	//sparse format the same as in deDoc.
		Option binsizeopt = Option.builder("binsize").argName("size").hasArg().required(true).desc("binsize of the input Hi-C data (kb)").build();
		Option outputfileopt = Option.builder("outputfile").argName("file").hasArg().desc("output Hi-C data file name").build();
		Option windowsizeopt = Option.builder("windowsize").argName("size").hasArg().desc("size of sliding window for deDoc2.s (bit), default 10(Mb)").build();
		Option maxTLDsizeopt = Option.builder("maxTLDsize").argName("size").hasArg().desc("maximum size of TLDs, to speed up deDoc2, default 10(Mb)").build();
		Option rwropt = new Option("rwr", "whether to do RWR imputation on input Hi-C, default false");
		Option rpopt = Option.builder("rp").argName("rp").hasArg().desc("restart probability of the RWR imputation process, default 0.5").build();
		Option NDIopt = new Option("NDI", "calculate the nromalized decoding information(NDI) of the input Hi-C data at current binsize");
		Option helpopt = new Option("help", "print usage information");
				
		options.addOption(inputfileopt);
		options.addOption(sparseformatopt);
		options.addOption(binsizeopt);
		options.addOption(outputfileopt);
		options.addOption(windowsizeopt);
		options.addOption(maxTLDsizeopt);
		options.addOption(rwropt);
		options.addOption(rpopt);
		options.addOption(NDIopt);
		options.addOption(helpopt);
		
		try {
			CommandLine line = parser.parse(options, args);
			if(line.hasOption("inputfile")){
				inputFileName = line.getOptionValue("inputfile");
			}
			if(line.hasOption("sparseformat")) {
				sparseformat = true;
			}
			if(line.hasOption("binsize")) {
				binsize = Integer.parseInt(line.getOptionValue("binsize"))*1000;
				if(binsize <= 0)
					throw new ParseException("the value of binsize should greater than 0.");
			}
			if(line.hasOption("outputfile")) {
				outputFileName = line.getOptionValue("outputfile");
			}
			else {
				outputFileName = inputFileName;
			}
			if(line.hasOption("windowsize")) {
				windowsize = new Double(Double.parseDouble(line.getOptionValue("windowsize"))*1000000).intValue();
				if(windowsize <= 0)
					throw new ParseException("the value of windowsize should greater than 0.");
			}
			if(line.hasOption("maxTLDsize")) {
				maxTLDsize = new Double(Double.parseDouble(line.getOptionValue("maxTLDsize"))*1000000).intValue();
				if(maxTLDsize <= 0)
					throw new ParseException("the value of maxTLDsize should greater than 0.");
			}
			if(line.hasOption("rwr")) {
				rwr = true;
			}
			if(line.hasOption("rp")) {
				rp = Double.parseDouble(line.getOptionValue("rp"));
				if(rp <= 0 || rp >= 1)
					throw new ParseException("the value of rp should be between 0 and 1 (not included).");
			}
			if(line.hasOption("NDI")) {
				NDI = true;
			}
			if(line.hasOption("help")) {
				HelpFormatter formatter = new HelpFormatter();
				formatter.printHelp( "deDoc2", options );
			}
			
			RealMatrix mtx;
			System.out.println("reading input data...");
			if(sparseformat) {
				mtx = Filerw.getMatrix_sparse(inputFileName);
			}
			else {
				mtx = Filerw.getMatrix(inputFileName);
			}
			RealMatrix mtx_org = mtx.copy();
			if(rwr) {
				System.out.println("performing rwr imputation...");
				mtx = Impute.rwr(mtx, rp, tol);
			}
			System.out.println("predicting TLDs...");
			DPAlgorithm algo = new DPAlgorithm(mtx, binsize, windowsize, maxTLDsize);
			algo.putTLD(outputFileName);
			if(NDI) {
				DPAlgorithm algo_NDI = new DPAlgorithm(mtx_org, binsize, windowsize, maxTLDsize);
				double SE1D = algo_NDI.calStrucEntro1D();
				double SE2D = algo_NDI.calStrucEntro2D(algo_NDI.getPartition2D());
				System.out.println("noralized decoding information of the input Hi-C matrix at binsize " + binsize + " is " + (SE1D-SE2D)/SE1D);
			}
			
		}catch(ParseException exp) {
			System.err.println( "Parsing failed.  Reason: " + exp.getMessage() );
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp( "deDoc2", options );
		}
	}
}
