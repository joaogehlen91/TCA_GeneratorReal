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
		double R = Integer.parseInt(algorithmParameters.get("R"));
		double d = Double.parseDouble(algorithmParameters.get("d"));
		double linkCapacities = Double.parseDouble(algorithmParameters.get("linkCapacities"));

		netPlan.removeAllNodes();
		netPlan.removeAllSRGs();
		
		R = fixNumberRegion(R);

		/* Generate node XY position table */
		for (int n = 0; n < N; n++)
		{
			double xCoord = xmin + (xmax - xmin) * r.nextDouble();
			double yCoord = ymin + (ymax - ymin) * r.nextDouble();
			netPlan.addNode(xCoord, yCoord, "Node " + n, null);
		}
		
		for (int i = 0; i < R; i++){
			
		}

		Set<Long> nodeIds = netPlan.getNodeIds();
		double dist_max = -Double.MAX_VALUE;
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

		return "Ok!";
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
	
	public Double capacityRegion(){
		
		
		return null;
	}
	
	public Double areaRegion(int R){
		
		
		return null;
	}
	
	public int fixNumberRegion(double r){
		boolean prime = true;
		for(int i=2; i<r; i++){
			if(r % i == 0){
				prime = false;
				break;
			}
		}

		if(prime){
			r = r + 1;
			return (int) r;
		}
		
		return (int) r;
	}
	
	
}