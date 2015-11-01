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
import com.net2plan.utils.Pair;
import com.net2plan.utils.Triple;


import java.util.*;

public class TCA_RealGenerator implements IAlgorithm
{	

	public static double lin, col;
	public static int p1 = 0, p2 = 0;
	List<Pair<Integer, Integer>> posicoes = new ArrayList<Pair<Integer, Integer>>();	
	List<List<Long>> idNodesRegioes = new ArrayList<List<Long>>();
	
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
		
		// função para arrumar o numero de regioes, caso seja primo.
		if(R > 2)
			R = fixNumberRegion(R);
	
		// x e y vão ser o tamanho do plano, ou seja o tamanho da matriz 
		int x = (int) (xmax-xmin);    
		int y = (int) (ymax-ymin);
		int plano[][] = new int[x][y];
		
		// zera todo o plano
		for (int[] row: plano) 		
			Arrays.fill(row, 0);
		
		/* 
		 * laço para encontrar p1 e p2. usado para separar as regioes de forma igual,
	     * conforme explicado no artigo do Professor. 
	     * p2 >= p1
	     * */
		for(int i = (R-1); i >= 1; i--){ 			
			if(((R % i) == 0) && isPrime(i)){	
				p1 = (R/i);
				p2 = i;
				break;
			}			
		}

		// as duas proximas linhas são para calcular a area de cadaReg para calcular a capacidade de cada regiao
		double areaReg = ((xmax-xmin)/p1) * ((ymax-ymin)/p2);
		int cap = (int)(areaReg/(d*d));
		
		// essas 10 linhas abaixo são para didivir melhor o numero de regioes, caso o plano nao seja um quadrado perfeito(x = y)
		// se o lado x for menor ou igual que y (x<=y) então plano vai ser dividido no lado do x por p1(que é o menor numero) e o
		// lado y por p2 já que é lado maior. Resumindo: divide o lado maior do plano com o numero maior de p(p1 ou p2).
		int sidex, sidey;      
		if((xmax-xmin) <= (ymax-ymin)){
			sidex = (int) ((xmax-xmin)/p1);
			sidey = (int) ((ymax-ymin)/p2);
			lin = p1; col = p2;
		}else{
			sidex = (int) ((xmax-xmin)/p2);
			sidey = (int) ((ymax-ymin)/p1);
			lin = p2; col = p1;
		}
		
		// matriz de R linhas para guardar onde inicia e termina cada região.
		// xi, yi, xf, yf
		//  0,  0,  2,  2
		//	0, 	2,  4,  4
		int indiceRegiao[][] = new int[R][4];
		int nreg = 0;
		for(int i = 0; i < lin; i++){
			for(int j = 0; j < col; j++){
				indiceRegiao[nreg][0] = sidex * i;
				indiceRegiao[nreg][1] = sidey * j;
				indiceRegiao[nreg][2] = sidex * (i+1);
				indiceRegiao[nreg][3] = sidey * (j+1);
				nreg++;
			}
		}
		
//		System.out.println("Lin = "+lin);
//		System.out.println("Col = "+col);
//		System.out.println("Capacidade de cada regiao = "+cap);
//		
//		for (int[] is : indiceRegiao) {
//			for (int i : is) {
//				System.out.print("\t"+i);
//			}
//			System.out.println();
//		}
		
		// laço que vai percorrer as regiões para posteriormente inserir os nodos.
		int nodosForReg = 0;
		idNodesRegioes.clear();
		for (int i = 0; i < R; i++) {
				
			// pega as posiçoes de inicio e fim da regiao i.
			int xi = indiceRegiao[i][0];		
			int yi = indiceRegiao[i][1];
			int xf = indiceRegiao[i][2];
			int yf = indiceRegiao[i][3];
			
			// preenche o vetor de posições com todas as posições da região i.
			posicoes.clear();			
			for (int ir = xi; ir < xf; ir++) {
				for (int jr = yi; jr < yf; jr++) {
					// * fazer um teste para nao pegar posição do plnao já ocupada
					posicoes.add(Pair.of(ir, jr));
				}
			}
			
			// sorteia quantos nodos vão ser inseridos na região R
			// 'i == 0' serve para evitar que muitos nos ficassem na primeira regiao.
			if(N == 1){
				nodosForReg = 1;
			}else{
				nodosForReg = r.nextInt(Math.min((i == 0 ? N/2 : N), (int)cap)); 
			}
			
			// laço para inserir cada nodo k na região i.
			List<Long> idNodesReg = new ArrayList<Long>();
			for(int k = 0; k < nodosForReg; k++){
				if(posicoes.size() > 0){
					int posicaoSorteada = r.nextInt(posicoes.size());
					int pi = posicoes.get(posicaoSorteada).getFirst();
					int pj = posicoes.get(posicaoSorteada).getSecond();
					
					if(plano[pi][pj] == 0){  // se está livre a posicao, mais um teste para garantir
						//netPlan.addNode(pi, pj, k+" "+i, null);
						idNodesReg.add(netPlan.addNode(pi, pj, k+" "+i, null));
						posicoes.remove(posicaoSorteada);
						plano = insereNodo(pi, pj, plano, x, y, (int) d);
						N--;
					}
				}
			}
			idNodesRegioes.add(idNodesReg);
			
			if(N==0) break;
		}
		
		// eu odeio java !!!!!
		List<List<Long>> idNodesRegioesAux = new ArrayList<List<Long>>();
		for (List<Long> list : idNodesRegioes) {
			List<Long> idNodesRegAux = new ArrayList<Long>(list);
			idNodesRegioesAux.add(idNodesRegAux);
		}
		
		/* daqui pra baixo é pra inserir os links(arestas)  */
		double dist_max = -Double.MAX_VALUE;
		for (List<Long> reg : idNodesRegioes) {
			
			//double dist_max = -Double.MAX_VALUE;
			for (long destinationNodeId : reg)
			{
				for (long originNodeId : reg)
				{
					if (originNodeId >= destinationNodeId) break;

					double dist = netPlan.getNodePairEuclideanDistance(originNodeId, destinationNodeId);
					if (dist > dist_max) dist_max = dist;
				}
			}
			
			
			int i = 0;
			List<Long> regOrdenada = new ArrayList<Long>();
			if(!reg.isEmpty()) regOrdenada.add(reg.get(i));
			//System.out.println(reg);
			while(reg.size()>1){
				int nexti = 0;
				long originNodeId = reg.get(i), idD = 0;
				double pmax = 0, dist = 0;
				for (int j = 0; j < reg.size(); j++) {
					long destinationNodeId = reg.get(j);
					if (originNodeId == destinationNodeId) continue;
					
					dist = netPlan.getNodePairEuclideanDistance(originNodeId, destinationNodeId);
					double p = alpha * Math.exp(-dist / (beta * dist_max));
					
					if(p > pmax) {
						pmax = p;
						idD = destinationNodeId;
						nexti = j;
					}
				}
				regOrdenada.add(idD);
				reg.remove(i);
				
				if(i > nexti){
					i = nexti;
				}else{
					i = (nexti -1);
				}
			}
			
			for (int j = 1; j < regOrdenada.size(); j++) {
				long origin = regOrdenada.get(j-1), destination = regOrdenada.get(j);
				if (origin == destination) continue;
				netPlan.addLink(origin, destination, linkCapacities, 10, null);
			}
			if(regOrdenada.size()>1)
				netPlan.addLink(regOrdenada.get(0), regOrdenada.get(regOrdenada.size()-1), linkCapacities, 10, null);
			
		}
		
		System.out.println(idNodesRegioesAux);
		
		
		//falta fazer essa parte: inserir links entre regioes
		for (List<Long> reg1 : idNodesRegioesAux) {
			long noOrigem = 0, noDestino = 0;
			double pmax = 0;
			for (Long no1 : reg1) {
				for (List<Long> reg2 : idNodesRegioesAux) {
					if(reg1 == reg2 ) continue;
					
					for (Long no2 : reg2) {
							
						double dist = netPlan.getNodePairEuclideanDistance(no1, no2);
						double p = alpha * Math.exp(-dist / (beta * dist_max));
						
						if(p > pmax) {
							pmax = p;
							noOrigem = no1;
							noDestino = no2;
						}
					}	
				}
			}
			// aqui add link
			if (noDestino != noOrigem) {
				netPlan.addLink(noOrigem, noDestino, linkCapacities, 10, null);
			}
			
		}

		
//		for (int l = 0; l < 100; l++) {
//			for (int c = 0; c < 100; c++) {
//				System.out.print(plano[l][c]);
//			}
//			System.out.println();
//		}
		
				
		return "Ok! | Nodos não inseridos: " + N;
	}

	private int[][] insereNodo(int pi, int pj, int[][] plano, int x, int y, int d) {
		plano[pi][pj] = 2;
		// seleciona a area para colocar a sombra
		int lini = Math.max(0, pi-d);
		int lfim = Math.min(x-1, pi+d);
		int cini = Math.max(0, pj-d);
		int cfim = Math.min(y-1, pj+d);
		
		// coloca a sombra e remove do vetor de posicoes
		for (int l = lini; l <= lfim; l++) {
			for (int c = cini; c <= cfim; c++) {
				posicoes.remove(Pair.of(l, c));
				if(l != pi || c != pj)
					plano[l][c] = 1;
			}
		}
		return plano;
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
		algorithmParameters.add(Triple.of("d", "6", "Minimum distance between two nodes"));
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
	
	
	public int fixNumberRegion(int r){
		if(isPrime(r))
			r = r + 1;
		return (int) r;
	}
	
	public boolean isPrime(int n){
		for(int i=2; i<n; i++){
			if(n % i == 0)
				return false;
		}
		return true;
	}

}