package tool;

import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.BlockRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;

public class Impute {
	public static RealMatrix sqrtVC(RealMatrix mtx) {
		mtx = mtx.copy();
		int dim = mtx.getRowDimension();
		RealMatrix rowSum = MatrixUtils.createRealMatrix(dim, 1);
		for(int i=0; i<dim; i++)
			rowSum.setEntry(i, 0, 1);
		rowSum = mtx.multiply(rowSum);
		RealMatrix diag = MatrixUtils.createRealDiagonalMatrix(rowSum.getColumn(0));
		for(int i=0; i<dim; i++) {
			double valuei = diag.getEntry(i, i);
			if(valuei == 0)
				diag.setEntry(i, i, 0);
			else
				diag.setEntry(i, i, 1/Math.sqrt(diag.getEntry(i, i)));
		}
//		RealMatrix sqrtVCMtx = diag.multiply(mtx.multiply(diag));
		RealMatrix sqrtVCMtx = mtx;
		for(int i=0; i<dim; i++) {
			for(int j=0; j<dim; j++) {
				sqrtVCMtx.setEntry(i, j, mtx.getEntry(i, j)*diag.getEntry(i, i)*diag.getEntry(j, j));
			}
		}
		mtx = null;
		sqrtVCMtx = symmetrizeMatrix(sqrtVCMtx);
		return sqrtVCMtx;
	}
	
	public static RealMatrix symmetrizeMatrix(RealMatrix mtx) {
		int dim = mtx.getRowDimension();
		for(int i=0;i<dim;i++)
			for(int j=i+1;j<dim;j++) {
				double weight = (mtx.getEntry(i, j) + mtx.getEntry(j, i)) / 2;
				mtx.setEntry(i, j, weight);
				mtx.setEntry(j, i, weight);
			}
		return mtx;
	}
	
	public static RealMatrix convo(RealMatrix mtx) {
		mtx = mtx.copy();
		double[][] filter = {{1,1,1},{1,2,1},{1,1,1}};
		int numOfSteps = 1;
		int dim = mtx.getRowDimension();
		RealMatrix convoMtx = new BlockRealMatrix(dim, dim);
		for(int i=0; i<dim; i++) {
			for(int j=0; j<dim; j++) {
				double valueij = 0;
				for(int m=-numOfSteps;m<=numOfSteps;m++) {
					for(int n=-numOfSteps;n<numOfSteps;n++) {
						if(i+m>=0 && i+m<dim && j+n>=0 && j+n<dim)
							valueij = valueij + filter[numOfSteps+m][numOfSteps+n] * mtx.getEntry(i+m, j+n);
					}
				}
				convoMtx.setEntry(i, j, valueij);
			}
		}
		return convoMtx;
	}
	
	public static RealMatrix convo(RealMatrix mtx, double[][] filter) {
		mtx = mtx.copy();
		int numOfSteps = filter.length / 2;
		if(filter[0].length != filter.length)
			throw new RuntimeException("filter is not square!");
		int dim = mtx.getRowDimension();
		RealMatrix convoMtx = new BlockRealMatrix(dim, dim);
		for(int i=0; i<dim; i++) {
			for(int j=0; j<dim; j++) {
				double valueij = 0;
				for(int m=-numOfSteps;m<=numOfSteps;m++) {
					for(int n=-numOfSteps;n<numOfSteps;n++) {
						if(i+m>=0 && i+m<dim && j+n>=0 && j+n<dim)
							valueij = valueij + filter[numOfSteps+m][numOfSteps+n] * mtx.getEntry(i+m, j+n);
					}
				}
				convoMtx.setEntry(i, j, valueij);
			}
		}
		convoMtx = symmetrizeMatrix(convoMtx);
		return convoMtx;
	}
	
	public static RealMatrix rwr(RealMatrix mtx, double alpha, double tol) {
		int dim = mtx.getRowDimension();
		RealMatrix W = sqrtVC(mtx);
		RealMatrix Q1 = MatrixUtils.createRealIdentityMatrix(dim);
		RealMatrix I = MatrixUtils.createRealIdentityMatrix(dim);
		RealMatrix Qnext;
		int iteration = 0;
		while(true) {
			Qnext = (Q1.multiply(W).scalarMultiply(1-alpha)).add(I.scalarMultiply(alpha));
			double norm2 = (Qnext.subtract(Q1)).getFrobeniusNorm();
			iteration ++;
//			System.out.println("iteration:\t" + iteration);
			Q1 = Qnext;
			if(norm2 <= tol) {
				System.out.println("iteration:\t" + iteration);
				break;
			}
			if(iteration >= 30) {
				System.err.println("random walk process didn't converge!");
				break;
			}
		}
		RealMatrix rwrMtx = symmetrizeMatrix(Q1);
		return rwrMtx;
	}
	
	public static RealMatrix rwr(RealMatrix mtx, double alpha) {
		int dim = mtx.getRowDimension();
		RealMatrix W = sqrtVC(mtx);
		RealMatrix Q1 = MatrixUtils.createRealIdentityMatrix(dim);
		RealMatrix I = MatrixUtils.createRealIdentityMatrix(dim);
		RealMatrix Qnext;
		int iteration = 0;
		while(true) {
			Qnext = (Q1.multiply(W).scalarMultiply(1-alpha)).add(I.scalarMultiply(alpha));
			double norm2 = (Qnext.subtract(Q1)).getFrobeniusNorm();
			iteration ++;
//			System.out.println("iteration:\t" + iteration);
			Q1 = Qnext;
			if(norm2 <= 0.01) {
				System.out.println("iteration:\t" + iteration);
				break;
			}
			if(iteration >= 50) {
				System.err.println("random walk process didn't converge!");
				break;
			}
		}
		RealMatrix rwrMtx = symmetrizeMatrix(Q1);
		return rwrMtx;
	}
	
	public static RealMatrix scHiCluster(RealMatrix mtx, double alpha) {
		RealMatrix convoMtx = convo(mtx);
		RealMatrix rwrMtx = rwr(convoMtx, alpha);
		return rwrMtx;
	}
}
