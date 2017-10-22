import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Scanner;

public class anchor_alignment {
	
	static int K;
	static int p1;
	static int p2;
	static int g;
	static int s;
	static String x;
	static String y;
	static int[] anchor_x;
	static int[] anchor_y;
	static ArrayList<Pair> anchorPair;
	static double[][] M;
	static double[][] Ix;
	static double[][] Iy;
	

	public static void main(String[] args) throws FileNotFoundException {
		// TODO Auto-generated method stub
		
		File sourceText = new File(args[0]);
		Scanner scnr = new Scanner(sourceText);
		anchorPair = new ArrayList<Pair>();
		
		K = scnr.nextInt(); //number of anchor pair
		p1 = scnr.nextInt(); //score for matched pairs
		p2 = scnr.nextInt(); //score for mismatched pairs
		g = scnr.nextInt(); //constant for gap penalty
		s = scnr.nextInt(); //linear coefficient for gap penalty
		scnr.nextLine();
		
		x = scnr.nextLine();
		y = scnr.nextLine();
		
		M = new double [x.length()+1][y.length()+1];
		Ix = new double [x.length()+1][y.length()+1];
		Iy = new double [x.length()+1][y.length()+1];
		
		if(K>0)
		{
		    anchor_x = new int[K];
		    anchor_y = new int[K];
    		for(int i = 0; i < K; i++)
    		{
    			anchor_x[i] = scnr.nextInt();
    			anchor_y[i] = scnr.nextInt();
    			anchorPair.add(new Pair(anchor_x[i], anchor_y[i]));
    			//System.out.println(anchor_x[i]+" "+anchor_y[i]);
	    	}
		
    		scnr.nextLine();
		}
		
//		System.out.println(K+" "+p1+" "+p2+" "+g+" "+s);
//		System.out.println(x);
//		System.out.println(y);
		
		fillMatrix();
		
//		System.out.println();
//		System.out.println("M");
//		for(int i = 0; i < M.length; i++)
//		{
//			for(int j = 0; j < M[0].length; j++)
//			{
//				System.out.print(M[i][j]+ " ");
//			}
//			System.out.println();
//		}
//		
//		System.out.println();
//		System.out.println("Ix");
//		for(int i = 0; i < Ix.length; i++)
//		{
//			for(int j = 0; j < Ix[0].length; j++)
//			{
//				System.out.print(Ix[i][j]+ " ");
//			}
//			System.out.println();
//		}
//		
//		System.out.println();
//		System.out.println("Iy");
//		for(int i = 0; i < Iy.length; i++)
//		{
//			for(int j = 0; j < Iy[0].length; j++)
//			{
//				System.out.print(Iy[i][j]+ " ");
//			}
//			System.out.println();
//		}
		Seq align = new Seq();
		if(K != 0)
		{
			align.concat(trace(0,0,anchorPair.get(0).x,anchorPair.get(0).y));
		
			for(int k = 0; k < K-1; k++)
			{
				align.concat(trace(anchorPair.get(k).x,anchorPair.get(k).y,
						anchorPair.get(k+1).x,anchorPair.get(k+1).y));
			}
			
			align.concat(trace(anchorPair.get(K-1).x,anchorPair.get(K-1).y,M.length-1,M[0].length-1));
			//System.out.println("("+anchorPair.get(K-1).x+","+anchorPair.get(K-1).y+") to ("+(M.length-1)+","+(M[0].length-1)+")");
		}
		else
		{
			align = trace(0,0,M.length-1,M[0].length-1);
		}
		
		System.out.print(align);
		
		

	}
	
	private static int S(int i, int j)
	{
		if(x.charAt(i-1) == y.charAt(j-1))
			return p1;
		else
			return (-p2);
	}
	
	private static int w(int k)
	{
		if(k == 0)
			return 0;
		
		return (-g-s*k);
	}
	
	
	private static void fillMatrix()
	{
		//initialize
		M[0][0] = (double) 0;
		Ix[0][0] = (double)w(0);
		Iy[0][0] = (double)w(0);
		for(int i = 1; i <= x.length(); i++)
		{
			M[i][0] = Double.NEGATIVE_INFINITY;
			Ix[i][0] = (double)w(i);
			Iy[i][0] = Double.NEGATIVE_INFINITY;
		}
		for(int j = 1; j <= y.length(); j++)
		{
			M[0][j] = Double.NEGATIVE_INFINITY;
			Ix[0][j] = Double.NEGATIVE_INFINITY;
			Iy[0][j] = (double)w(j);
		}
		
		//fill in rest of matrices
		for(int i = 1; i < M.length; i++)
		{
			for(int j = 1; j < M[0].length; j++)
			{
				for(int k = 0; k < anchorPair.size(); k++)
				{
					if(anchorPair.get(k).match(i, j))
					{
						M[i][j] = Math.max(Math.max(M[i-1][j-1]+S(i,j), Ix[i-1][j-1]+S(i,j)),Iy[i-1][j-1]+S(i,j));
						Ix[i][j] = M[i][j];
						Iy[i][j] = M[i][j];
					}
				}
				
				//System.out.println(i+" "+j);
				M[i][j] = Math.max(Math.max(M[i-1][j-1]+S(i,j), Ix[i-1][j-1]+S(i,j)),Iy[i-1][j-1]+S(i,j));
				Ix[i][j] = Math.max(M[i-1][j]+w(1), Ix[i-1][j]-s);
				Iy[i][j] = Math.max(M[i][j-1]+w(1), Iy[i][j-1]-s);
			}
		}
		
		
	}
	
	private static Seq trace(int x0, int y0, int i, int j){
		
		char tracker;
		if(M[i][j] == max3(i,j))
		{
			tracker = 'm';
		}
		else if (Ix[i][j] == max3(i,j))
		{
			tracker = 'x';
		}
		else
		{
			tracker = 'y';
		}
		Seq align = new Seq();
		
		while(i >= x0 && j >= y0)
		{
			if(i == x0 && j == y0)
				break;
			if(i == x0)
				tracker = 'y';
			if(j == y0)
				tracker = 'x';
			
			switch(tracker){
			case 'm':
				if(M[i-1][j-1] == max3(i-1,j-1))
				{
					tracker = 'm';
				}
				else if (Ix[i-1][j-1] == max3(i-1,j-1))
				{
					tracker = 'x';
				}
				else
				{
					tracker = 'y';
				}
				//System.out.println("M: "+i+" "+j+" " + tracker);
				align.add(x.charAt(i-1),y.charAt(j-1));  //POSSIBLE BUG
				i = i-1;
				j = j-1;
				break;
				
			case 'x':
				if(Ix[i][j] == M[i-1][j]+w(1))
				{
					tracker = 'm';
				}
				else
				{
					tracker = 'x';
				}
				//System.out.println("X gap: "+i+" "+j+" " + tracker);
				align.add(x.charAt(i-1), '_');  //POSSIBLE BUG
				i = i-1;
				break;
				
			case 'y':
				if(Iy[i][j] == M[i][j-1]+w(1))
				{
					tracker = 'm';
				}
				else
				{
					//System.out.println(Iy[i][j]+" +"+s+" = "+Iy[i][j-1]+"; w(1)="+w(1));
					tracker = 'y';
				}
				//System.out.println("Y gap: "+i+" "+j+" " + tracker);
				align.add('_', y.charAt(j-1));  //POSSIBLE BUG
				j = j-1;
				break;
			}
			
		}
		return align;
	}
	
	private static double max3(int i, int j) {
		return Math.max(Math.max(M[i][j], Ix[i][j]), Iy[i][j]);
	}
	
	
}

class Node {
	char t;
	double v;
	
	Node(double value, char tracker)
	{
		this.v = value;
		this.t = tracker;
	}
}

class Seq {
	String x_seq;
	String y_seq;
	
	
	public Seq()
	{
		x_seq = "";
		y_seq = "";
	}
	
	public void add(char x, char y)
	{
		x_seq = x + x_seq;
		y_seq = y + y_seq;
	}
	
	public void concat(Seq seq)
	{
		this.x_seq = this.x_seq + seq.x_seq;
		this.y_seq = this.y_seq + seq.y_seq;
	}
	
	public String toString(){
		return x_seq+"\n"+y_seq+"\n";
	}
}


class Pair {
	
	int x;
	int y;
	
	public Pair(int x, int y){
		this.x = x;
		this.y = y;
	}
	
	public boolean match(int i, int j)
	{
		if(i == x && j == y)
			return true;
		else
			return false;
	}
}

	
