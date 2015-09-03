/*******************************************************************************
 * Copyright (c) 2013-2015 Pablo Pavon-Marino, Jose-Luis Izquierdo-Zaragoza.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Lesser Public License v3
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/lgpl.html
 *
 * Contributors:
 *     Pablo Pavon-Marino, Jose-Luis Izquierdo-Zaragoza - initial API and implementation
 ******************************************************************************/
package com.net2plan.examples.netDesignAlgorithm.tca;

import com.net2plan.interfaces.networkDesign.IAlgorithm;
import com.net2plan.interfaces.networkDesign.NetPlan;
import com.net2plan.utils.Triple;
import java.util.*;

public class TCA_RealGenerator implements IAlgorithm
{	

	public static Map<String, String> region = new HashMap<String, String>();
	public static double lin, col;
	
	@Override
	public String executeAlgorithm(NetPlan netPlan, Map<String, String> algorithmParameters, Map<String, String> net2planParameters)
	{
		Random r = new Random();
		

		int N = Integer.parseInt(algorithmParameters.get("N"));
		double alpha = Double.parseDouble(algorithmParameters.get("alpha"));
		double beta = Double.parseDouble(algorithmParameters.get("beta"));
		double xmax = Double.parseDouble(algorithmParameters.get("xmax"));
		double xmin = Double.parseDouble(algorithmParameters.get("xmin"));
		double ymax = Double.parseDouble(algorithmParameters.get("ymax"));
		double ymin = Double.parseDouble(algorithmParameters.get("ymin"));
		int R = Integer.parseInt(algorithmParameters.get("R"));
		double d = Double.parseDouble(algorithmParameters.get("d"));
		double linkCapacities = Double.parseDouble(algorithmParameters.get("linkCapacities"));

		netPlan.removeAllNodes();
		netPlan.removeAllSRGs();
		netPlan.reset();
		
		R = fixNumberRegion(R);
		int cap = capacityRegion(R, xmin, xmax, ymin, ymax, d);
		double sidex = Double.parseDouble(region.get("sidex"));
		double sidey = Double.parseDouble(region.get("sidey"));

		int ninreg = 0;
				
		/* laço para percorrer as regioes e inserir os nodos em cada regiao
		 * 'sidex' e 'sidey' são os tamanhos dos lados de cada regiao(para poder navegar no plano);
		 * 'ninreg' guarda a quantidade de nos que vão são inseridos na região;
		*/
		int reg = 0; 		
		while(N > 1){
			for(int x = 0; x < lin; x++){
				for(int y = 0; y < col; y++){
					ninreg = r.nextInt(Math.min(N, (int)cap));
					for(int k = 0; k < ninreg; k++){
						double xCoord2 = randomRegion(sidex*x, sidex*(x+1));
						double yCoord2 = randomRegion(sidey*y, sidey*(y+1));
						netPlan.addNode(xCoord2, yCoord2, "Node " + k +"reg "+reg, null);
					}
					N = (int)(N - ninreg);
					reg++;
				}
			}
		}
		
		// se ainda tem um nó para inserir, insere na ultima região, isso ocorre por causa do sorteio, pode ser que o ultimo não sortei.
		if(N > 0){
			double xCoord2 = randomRegion(sidex*(lin-1), sidex*lin);
			double yCoord2 = randomRegion(sidey*(col-1), sidey*col);
			netPlan.addNode(xCoord2, yCoord2, "Node extra" +"reg "+sidex*(lin-1)+" "+sidey*col+" "+sidex+" "+sidey, null);
		}
		

		Set<Long> nodeIds = netPlan.getNodeIds();
		//double dist_max = -Double.MAX_VALUE;
		double dist_max = d;
		for (long destinationNodeId : nodeIds)
		{
			for (long originNodeId : nodeIds)
			{
				if (originNodeId >= destinationNodeId) break;

				double dist = netPlan.getNodePairEuclideanDistance(originNodeId, destinationNodeId);
				if (dist > dist_max) dist_max = dist;
			}
		}

		/* Generate a directed link between each node pair with probability p = alpha * exp(-distance/(beta * max_distance)) */
		for (long destinationNodeId : nodeIds)
		{
			for (long originNodeId : nodeIds)
			{
				if (originNodeId >= destinationNodeId) break;
				double dist = netPlan.getNodePairEuclideanDistance(originNodeId, destinationNodeId);
				double p = alpha * Math.exp(-dist / (beta * dist_max));

				if (r.nextDouble() < p)
					netPlan.addLink(originNodeId, destinationNodeId, linkCapacities, dist, null);
			}
		}
		
		return "Ok!" + "  " + region.get("area")+"N "+N;
	}

	@Override
	public List<Triple<String, String, String>> getParameters()
	{
		List<Triple<String, String, String>> algorithmParameters = new LinkedList<Triple<String, String, String>>();
		algorithmParameters.add(Triple.of("N", "30", "Number of nodes"));
		algorithmParameters.add(Triple.of("alpha", "0.4", "'alpha' factor"));
		algorithmParameters.add(Triple.of("beta", "0.4", "'beta' factor"));
		algorithmParameters.add(Triple.of("xmax", "100", "Right endpoint for the x-axis"));
		algorithmParameters.add(Triple.of("xmin", "0", "Left endpoint for the x-axis"));
		algorithmParameters.add(Triple.of("ymax", "100", "Upper endpoint for the y-axis"));
		algorithmParameters.add(Triple.of("ymin", "0", "Lower endpoint for the y-axis"));
		algorithmParameters.add(Triple.of("R", "6", "Regions"));
		algorithmParameters.add(Triple.of("d", "2", "Minimum distance between two nodes"));
		algorithmParameters.add(Triple.of("linkCapacities", "100", "The capacities to set in the links"));
		return algorithmParameters;
	}

	@Override
	public String getDescription()
	{
		String NEWLINE = String.format("%n");

		StringBuilder aux = new StringBuilder();
		aux.append("Descricao teste");
		//aux.append("This algorithm implements the random network topology generator introduced in Waxman (1988). The Waxman's generator is a geographic model for the growth of a network. In this model nodes are uniformly distributed in a given area and links are added according to probabilities that depend on the distances between the nodes. The probability to have a link between nodes i and j is given by:");
		aux.append(NEWLINE);
		aux.append(NEWLINE);
		//aux.append("P(i, j) = alpha * exp(-d/(beta * d_max)");
		aux.append(NEWLINE);
		aux.append(NEWLINE);
		//aux.append("where 0<alpha, beta<=1, d is the distance from i to j, and d_max is the maximum distance between any two nodes. An increase in the parameter alpha increases the probability of edges between any nodes in the network, while an increase in beta yields a larger ratio of long links to short links.");

		return aux.toString();
	}
	
	public double randomRegion(double min, double max){
		Random r = new Random();
		double coord = min + (max - min) * r.nextDouble();
		return coord;
	}
	
	// retorno a capacidade de cada região
	public Integer capacityRegion(int r, double xmin, double xmax, double ymin, double ymax, double d){
		double areaReg = areaRegion(r, xmin, xmax, ymin, ymax);
		return (int)(areaReg/(d*d));
	}
	
	public Double areaRegion(int r, double xmin, double xmax, double ymin, double ymax){
		double p1=0, p2=0;
		for(int i=(r-1); i>1; i--){
			if(((r % i) == 0) && isPrime(i)){
				p1 = (r/i);
				p2 = i;
				break;
			}			
		}
		
		double area = ((xmax-xmin)/p1) * ((ymax-ymin)/p2);
		
		region.put("area", String.valueOf(area));
		if((xmax-xmin) <= (ymax-ymin)){
			region.put("sidex", String.valueOf((xmax-xmin)/p1));
			region.put("sidey", String.valueOf((ymax-ymin)/p2));
			lin = p1; col = p2;
		}else{
			region.put("sidex", String.valueOf((xmax-xmin)/p2));
			region.put("sidey", String.valueOf((ymax-ymin)/p1));
			lin = p2; col = p1;
		}
		return area;
	}
	
	
	public int fixNumberRegion(int r){
		if(isPrime(r)){
			r = r + 1;
		}
		return (int) r;
	}
	
	public boolean isPrime(int n){
		for(int i=2; i<n; i++){
			if(n % i == 0){
				return false;
			}
		}
		return true;
	}
	
	
}
