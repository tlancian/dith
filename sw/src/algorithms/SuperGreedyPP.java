
package algorithms;

//import org.json.simple.JSONObject;
import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import java.util.Random;
import tools.FastAndSimpleIndexMinPQ;

/**
 * 
 * 
 */
public class SuperGreedyPP implements Serializable {

	private static final long serialVersionUID = 7802057891779118703L;

	public static void main(String[] args) throws Exception {
		//
		// input: [method, G, A, B, lambda_1, lambda_2, file_prox, file_dist]

//		serializeGraphOnFile("/Users/ikki/DATASETS/NegDSG/DATASETS/Amazon/com_amazon_ungraph__FINAL__FIRST_10000_ROWS.tsv");
//		if (args.length >= 0) {
//			return;
//		}
//		exp_ON_MODEL_for_HARD_DOWN_IN_THE_HOLLOW("");
		System.out.println();
		System.out.println(Arrays.toString(args));
		System.out.println();
		//
//		args = new String[8];
//		args[0] = "/Users/ikki/Downloads/ISI/peradri"; 
//		args[1] = "twitch_gamers__4US.tsv";
//		args[2] = "empty_file.txt";
//		args[3] = "empty_file.txt";
//		args[4] = "0";
//		args[5] = "0";
//		args[6] = "10000";
//		args[7] = "0.99";
		//
//		exp_RANDOM_SEARCH_ON_BLADE(args[0], args[1], args[2], args[3]);
//		exp_BASELINES_ON_BLADE(args[0], args[1], args[2], args[3]);
//		exp_GRID_SEARCH_ON_BLADE(args[0], args[1], args[2], args[3]);
		//
		//
		//
		//
//		exp_GRID_SEARCH_ON_BLADE_with_MAX_lambdas(args[0], args[1], args[2], args[3]);
//		if (args.length >= 0) {
//			return;
//		}
		//
//		exp_ON_BLADE_with_FIXED_lambdas__4_CONVERGENCE(args[0], args[1], args[2], args[3], Double.parseDouble(args[4]),
//				Double.parseDouble(args[5]), Integer.parseInt(args[6]), Double.parseDouble(args[7]));
//		if (args.length >= 0) {
//			return;
//		}
		//
//		args = new String[13];
//		args[0] = "DITH";
////		args[1] = "/Users/ikki/Downloads/ISI/peradri";
//		args[1] = "/Users/ikki/Downloads/QUALITY__OUTPUT_FINAL";
////		args[2] = "ego_facebook__4US.tsv";
//		args[2] = "leadersdebate__4US.tsv";
//		args[3] = "[0, 1 , 2]";
//		args[4] = "[0, 1 , 2]"; 
//		args[5] = "1";
//		args[6] = "1";
//		args[7] = "/Users/ikki/Downloads/ISI/peradri/empty_file.txt";
//		args[8] = "/Users/ikki/Downloads/ISI/peradri/empty_file.txt";
//		args[9] = "10000";
//		args[10] = "0.99";
//		args[11] = "/Users/ikki/Downloads/ISI/peradri/output";
//		args[12] = "42";
		String method = args[0].toUpperCase();
		String common_directory = args[1];
		String graph_file_name = args[2];
		String A = args[3];
		String B = args[4];
		double lambda_1 = Double.parseDouble(args[5]);
		double lambda_2 = Double.parseDouble(args[6]);
		String input_file_name__PROX = args[7];
		String input_file_name__DIST = args[8];
		int T = Integer.parseInt(args[9]);
		double gamma_approximation_factor = Double.parseDouble(args[10]);
		String output_directory = args[11];
		String exp_instance_id = args[12];
		if (method.equalsIgnoreCase("DITH")) {
			exp_ON_BLADE_with_FIXED_lambdas__4_QUALITY(common_directory, graph_file_name, input_file_name__PROX,
					input_file_name__DIST, lambda_1, lambda_2, T, gamma_approximation_factor, A, B, output_directory,
					exp_instance_id, "DITH");
		} else if (method.equalsIgnoreCase("DITH-1")) {
			T = 1;
			exp_ON_BLADE_with_FIXED_lambdas__4_QUALITY(common_directory, graph_file_name, input_file_name__PROX,
					input_file_name__DIST, lambda_1, lambda_2, T, gamma_approximation_factor, A, B, output_directory,
					exp_instance_id, "DITH-1");
		} else {
			exp_ON_BLADE_with_FIXED_lambdas__4_QUALITY_BASELINE(common_directory, graph_file_name,
					input_file_name__PROX, input_file_name__DIST, lambda_1, lambda_2, T, gamma_approximation_factor, A,
					B, output_directory, exp_instance_id, method);
		}
		if (args.length >= 0) {
			return;
		}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
		//
		int[] indexes = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600,
				700, 800, 1000 };
		for (Integer iiiindex : indexes) {

			String[] new_args = new String[4];
			new_args[0] = "/Users/ikki/Downloads/ANECDOTAL_EVIDENCE/greek_parliament";
			new_args[1] = "greek_parliament__NO_ZERO_WEIGHT_ON_EDGES.tsv";
			new_args[2] = "greek_parliament_query_node_1__MSSP_01.tsv";
			new_args[3] = "greek_parliament_query_node_2__1_minus_MSSP_01.tsv";
			//
			new_args[0] = "/Users/ikki/Dropbox/NegDSG/sw/MODEL_3/SCALABILITY/" + iiiindex;
			new_args[1] = "sbm.tsv";
			new_args[2] = "T__1_minus_MSSP_01.tsv";
			new_args[3] = "D__MSSP_01.tsv";
			args = new_args;

//		exp_BASELINES_ON_BLADE(args[0], args[1], args[2], args[3]);
			// exp_GRID_SEARCH_ON_BLADE_with_MAX_lambdas(args[0], args[1], args[2],
			// args[3]);
			exp_SUPER_GREEDY_PP_ON_BLADE_with_lambdas_1_1_1(args[0], args[1], args[2], args[3]);
		}

//		exp_ON_TWEETTER("");
//		exp_ON_MODEL("");
//		exp_12("flowless");
//		exp_11("flowless");
//		exp_10("flowless");
		if (args.length >= 0) {
			return;
		}

		return;
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	/**
	 * 
	 * @param c_commoon_path
	 * @param input_file_name
	 * @param all_c_input_file_name__set_T
	 * @param all_c_input_file_name__set_D
	 * @throws Exception
	 */
	protected static void exp_BASELINES_ON_BLADE(String c_commoon_path, String input_file_name,
			String all_c_input_file_name__set_T, String all_c_input_file_name__set_D) throws Exception {

		int prefix_index = c_commoon_path.lastIndexOf("/");
		String prefix = c_commoon_path.substring(prefix_index + 1);

		String c_input_file_name = c_commoon_path + "/" + input_file_name;
		BufferedWriter bw = new BufferedWriter(new FileWriter(
				c_commoon_path + "/BASELINES__" + prefix + "__" + input_file_name + "__" + all_c_input_file_name__set_T
						+ "__" + all_c_input_file_name__set_D + "__" + System.currentTimeMillis() + ".tsv"));
		//
		SuperGreedyPP neg_dsg_c = new SuperGreedyPP(c_input_file_name);
		neg_dsg_c.createGraph();
		neg_dsg_c.createSetsOfNodesToConsider(c_commoon_path + "/" + all_c_input_file_name__set_T,
				c_commoon_path + "/" + all_c_input_file_name__set_D, 0., 0.);
		//
		neg_dsg_c.lambda_0 = 1.;
		neg_dsg_c.lambda_1 = 0.;
		neg_dsg_c.lambda_2 = 0.;
		//
		String res_baseline_0 = neg_dsg_c.compute_baseline_MAX_DENSITY();// ok
		String res_baseline_1 = neg_dsg_c.compute_baseline_MAX_DISTANCE();// ok
		String res_baseline_2 = neg_dsg_c.compute_baseline_MAX_PROXIMITY();// ok
		String res_baseline_3 = neg_dsg_c.compute_baseline_MAX_DISTANCE_THEN_DENSITY();// ok
		String res_baseline_4 = neg_dsg_c.compute_baseline_MAX_PROXIMITY_THEN_DENSITY();// ok
		String res_baseline_5 = neg_dsg_c.compute_baseline_MAX_PROXIMITY_AND_DISTANCE_THEN_DENSITY();// ok
		//
		// bw.write(result);
		bw.write("res_baseline_MAX_DENSITY" + "\t" + res_baseline_0);
		bw.write("\n");
		bw.write("res_baseline_MAX_DISTANCE" + "\t" + res_baseline_1);
		bw.write("\n");
		bw.write("res_baseline_MAX_PROXIMIITY" + "\t" + res_baseline_2);
		bw.write("\n");
		//
		bw.write("baseline_MAX_DISTANCE_THEN_DENSITY" + "\t" + res_baseline_3);
		bw.write("\n");
		bw.write("baseline_MAX_PROXIMITY_THEN_DENSITY" + "\t" + res_baseline_4);
		bw.write("\n");
		//
		bw.write("baseline_MAX_PROXIMITY_AND_DISTANCE_THEN_DENSITY" + "\t" + res_baseline_5);
		bw.write("\n");
		bw.flush();

		//

		bw.close();
		return;
	}

	/////////////////////////////////////////////////////////////////////

	protected static void exp_ON_BLADE_with_FIXED_lambdas__4_CONVERGENCE(String c_commoon_path, String input_file_name,
			String all_c_input_file_name__set_T, String all_c_input_file_name__set_D, double l_1, double l_2, int T,
			double gamma_approximation_factor) throws Exception {

		int prefix_index = c_commoon_path.lastIndexOf("/");
		String prefix = c_commoon_path.substring(prefix_index + 1);

		String c_input_file_name = c_commoon_path + "/" + input_file_name;
		BufferedWriter bw = new BufferedWriter(new FileWriter(
				c_commoon_path + "/" + prefix + "__" + input_file_name + "__" + all_c_input_file_name__set_T + "__"
						+ all_c_input_file_name__set_D + "__" + System.currentTimeMillis() + ".tsv"));

		double lambda_1 = l_1;
		double lambda_2 = l_2;
//		int T = 10000;
		double min_appx_factor = gamma_approximation_factor;
//		double min_appx_factor = 1.;

		String result = "";

		Set<String> all_solutions_found_so_far = new HashSet<String>();

		SuperGreedyPP neg_dsg_c = new SuperGreedyPP(c_input_file_name);
		neg_dsg_c.createGraph();
		neg_dsg_c.createSetsOfNodesToConsider(c_commoon_path + "/" + all_c_input_file_name__set_T,
				c_commoon_path + "/" + all_c_input_file_name__set_D, lambda_1, lambda_2);
		//
		neg_dsg_c.lambda_0 = 1.;
		neg_dsg_c.lambda_1 = lambda_1;
		neg_dsg_c.lambda_2 = lambda_2;
		//
		System.out.println(neg_dsg_c.lambda_0 + "\t" + neg_dsg_c.lambda_1 + "\t" + neg_dsg_c.lambda_2);
		//
		System.out.println("-_-");
		result = neg_dsg_c.runSuperGreedyPlusPlusOnlySolutionValue_2(T, min_appx_factor, true);
		System.out.println("_-_");
		int last_index = result.lastIndexOf("\t");
		String c_solution = result.substring(last_index + 1);
		if (all_solutions_found_so_far.add(c_solution)) {
			System.out.println("Numer of all_solutions_found_so_far: " + all_solutions_found_so_far.size());
			bw.write(result);
			bw.write("\n");
			bw.flush();
		}

		bw.close();
		return;
	}

	protected static void exp_ON_BLADE_with_FIXED_lambdas__4_QUALITY(String c_commoon_path, String input_file_name,
			String all_c_input_file_name__set_T, String all_c_input_file_name__set_D, double l_1, double l_2, int T,
			double gamma_approximation_factor, String A, String B, String output_directory, String exp_instance_id,
			String method) throws Exception {

//		int prefix_index = c_commoon_path.lastIndexOf("/");
//		String prefix = c_commoon_path.substring(prefix_index + 1);

		String c_input_file_name = c_commoon_path + "/" + input_file_name;
		BufferedWriter bw = new BufferedWriter(new FileWriter(output_directory + "/" + "QUALITY_EXP_OUTPUT__"
				+ exp_instance_id + "__" + input_file_name + "__" + System.currentTimeMillis() + ".tsv"));

		double lambda_1 = l_1;
		double lambda_2 = l_2;
//		int T = 10000;
		double min_appx_factor = gamma_approximation_factor;
//		double min_appx_factor = 1.;

		String result = "";

		Set<String> all_solutions_found_so_far = new HashSet<String>();

		SuperGreedyPP algo = new SuperGreedyPP(c_input_file_name);
		algo.createGraph();
		algo.createSetsOfNodesToConsider(all_c_input_file_name__set_T, all_c_input_file_name__set_D, lambda_1,
				lambda_2);
		//
		algo.lambda_0 = 1.;
		algo.lambda_1 = lambda_1;
		algo.lambda_2 = lambda_2;
		//
		System.out.println(algo.lambda_0 + "\t" + algo.lambda_1 + "\t" + algo.lambda_2);
		//
		System.out.println("-_-");
		result = algo.runSuperGreedyPlusPlusOnlySolutionValue_2(T, min_appx_factor);
		System.out.println("_-_");
		int last_index = result.lastIndexOf("\t");
		String c_solution = result.substring(last_index + 1);
		if (all_solutions_found_so_far.add(c_solution)) {
			System.out.println("Numer of all_solutions_found_so_far: " + all_solutions_found_so_far.size());
			//
//			this.map__field__value.put("AvgDeg(Sol)", "" + density);
//			this.map__field__value.put("AvgProx(Sol)", "" + proximity_to_tau);
//			this.map__field__value.put("AvgDist(Sol)", "" + distance_from_danger);
//			this.map__field__value.put("OF(Sol)", "" + best_solution_so_far_value);
//			this.map__field__value.put("UB", "" + UB);
//			this.map__field__value.put("AppxFactor", "" + solution_quality);
//			this.map__field__value.put("t", "" + t);
//			this.map__field__value.put("ExecTimeMilliSec", "" + "" + (t_1 - t_0));
//			this.map__field__value.put("T", "" + T);
//			this.map__field__value.put("gamma", "" + min_appx_factor);
//			this.map__field__value.put("|Sol|", "" + best_solution_so_far_size);
//			this.map__field__value.put("Sol", subgraph_as_string);
			//
			result = "";
			result += method + "\t";
			result += input_file_name + "\t";
			result += A + "\t";
			result += B + "\t";
			result += algo.lambda_1 + "\t";
			result += algo.lambda_2 + "\t";
			result += algo.map__field__value.get("AvgDeg(Sol)") + "\t";
			result += algo.map__field__value.get("AvgProx(Sol)") + "\t";
			result += algo.map__field__value.get("AvgDist(Sol)") + "\t";
			result += algo.map__field__value.get("OF(Sol)") + "\t";
			result += algo.map__field__value.get("UB") + "\t";
			result += algo.map__field__value.get("AppxFactor") + "\t";
			result += algo.map__field__value.get("t") + "\t";
			result += algo.map__field__value.get("ExecTimeMilliSec") + "\t";
			result += algo.map__field__value.get("T") + "\t";
			result += algo.map__field__value.get("gamma") + "\t";
			result += algo.map__field__value.get("|Sol|") + "\t";
			result += algo.map__field__value.get("Sol");
			// output: [method, Graph_name, A, B, lambda_1, lambda_2, AvgDeg(Sol),
			// AvgProx(Sol), AvgDist(Sol), OF(Sol), UB, AppxFactor, t, ExecTimeMilliSec, T,
			// gamma, |Sol|, Sol]
			//
			System.out.println(result);
			bw.write(result);
			bw.write("\n");
			bw.flush();
		}

		bw.close();
		return;
	}

	protected static void exp_ON_BLADE_with_FIXED_lambdas__4_QUALITY_BASELINE(String c_commoon_path,
			String input_file_name, String all_c_input_file_name__set_T, String all_c_input_file_name__set_D,
			double l_1, double l_2, int T, double gamma_approximation_factor, String A, String B,
			String output_directory, String exp_instance_id, String baseline) throws Exception {

//		int prefix_index = c_commoon_path.lastIndexOf("/");
//		String prefix = c_commoon_path.substring(prefix_index + 1);

		String c_input_file_name = c_commoon_path + "/" + input_file_name;
		BufferedWriter bw = new BufferedWriter(new FileWriter(output_directory + "/" + "QUALITY_EXP_OUTPUT__"
				+ exp_instance_id + "__" + input_file_name + "__" + System.currentTimeMillis() + ".tsv"));

		double lambda_1 = l_1;
		double lambda_2 = l_2;
//		int T = 10000;
//		double min_appx_factor = gamma_approximation_factor;
//		double min_appx_factor = 1.;

		String result = "";

		Set<String> all_solutions_found_so_far = new HashSet<String>();

		SuperGreedyPP algo = new SuperGreedyPP(c_input_file_name);
		algo.createGraph();
		algo.createSetsOfNodesToConsider(all_c_input_file_name__set_T, all_c_input_file_name__set_D, lambda_1,
				lambda_2);
		//
		algo.lambda_0 = 1.;
		algo.lambda_1 = lambda_1;
		algo.lambda_2 = lambda_2;
		//
		System.out.println(algo.lambda_0 + "\t" + algo.lambda_1 + "\t" + algo.lambda_2);
		//
		System.out.println("-_-");
		long t0 = System.currentTimeMillis();
		if (baseline.equalsIgnoreCase("DS")) {
			result = algo.compute_baseline_MAX_DENSITY();
		} else if (baseline.equalsIgnoreCase("Ego-Prox-DS")) {
			result = algo.compute_baseline_MAX_PROXIMITY_THEN_DENSITY();
		} else if (baseline.equalsIgnoreCase("Hybrid-DS")) {
			result = algo.compute_baseline_MAX_PROXIMITY_AND_DISTANCE_THEN_DENSITY();
		} else {
			bw.close();
			return;
		}
		long t1 = System.currentTimeMillis();
		algo.map__field__value.put("ExecTimeMilliSec", "" + (t1 - t0));
		//
//		result = algo.runSuperGreedyPlusPlusOnlySolutionValue_2(T, min_appx_factor);
		System.out.println("_-_");
		int last_index = result.lastIndexOf("\t");
		String c_solution = result.substring(last_index + 1);
		if (all_solutions_found_so_far.add(c_solution)) {
			System.out.println("Numer of all_solutions_found_so_far: " + all_solutions_found_so_far.size());
			result = "";
			result += baseline + "\t";
			result += input_file_name + "\t";
			result += A + "\t";
			result += B + "\t";
			result += algo.lambda_1 + "\t";
			result += algo.lambda_2 + "\t";
			result += algo.map__field__value.get("AvgDeg(Sol)") + "\t";
			result += algo.map__field__value.get("AvgProx(Sol)") + "\t";
			result += algo.map__field__value.get("AvgDist(Sol)") + "\t";
			result += algo.map__field__value.get("OF(Sol)") + "\t";
			result += algo.map__field__value.get("UB") + "\t";
			result += algo.map__field__value.get("AppxFactor") + "\t";
			result += algo.map__field__value.get("t") + "\t";
			result += algo.map__field__value.get("ExecTimeMilliSec") + "\t";
			result += algo.map__field__value.get("T") + "\t";
			result += algo.map__field__value.get("gamma") + "\t";
			result += algo.map__field__value.get("|Sol|") + "\t";
			result += algo.map__field__value.get("Sol");
			// output: [method, Graph_name, A, B, lambda_1, lambda_2, AvgDeg(Sol),
			// AvgProx(Sol), AvgDist(Sol), OF(Sol), UB, AppxFactor, t, ExecTimeMilliSec, T,
			// gamma, |Sol|, Sol]
			//
			System.out.println(result);
			bw.write(result);
			bw.write("\n");
			bw.flush();
		}

		bw.close();
		return;
	}

	protected static double[] get_max_ALL_lambdas(String c_commoon_path, String input_file_name,
			String all_c_input_file_name__set_T, String all_c_input_file_name__set_D, int number_of_distinct_lambdas)
			throws Exception {
		//
		double[] all_lambdas = new double[number_of_distinct_lambdas + 1];
		//
		String c_input_file_name = c_commoon_path + "/" + input_file_name;
		//
		SuperGreedyPP dsg = new SuperGreedyPP(c_input_file_name);
		dsg.createGraph();
		dsg.createSetsOfNodesToConsider(c_commoon_path + "/" + all_c_input_file_name__set_T,
				c_commoon_path + "/" + all_c_input_file_name__set_D, 0, 0);
		//
		dsg.compute_baseline_MAX_DENSITY();
		double dsg_density = dsg.computeSubgraphDensity(dsg.output_best_V);
		DSG_density = dsg_density;
		//
		double MAX_lambda = 1. * dsg_density / 0.1;
		//
		double single_lambda_step = MAX_lambda / ((double) number_of_distinct_lambdas);
		for (int c_factor = 0; c_factor <= number_of_distinct_lambdas; c_factor++) {
			all_lambdas[c_factor] = c_factor * single_lambda_step;
		}
		//
		return all_lambdas;
	}

	protected static void exp_GRID_SEARCH_ON_BLADE_with_MAX_lambdas(String c_commoon_path, String input_file_name,
			String all_c_input_file_name__set_T, String all_c_input_file_name__set_D) throws Exception {

//		int number_of_distinct_lambdas = 50;
//		int number_of_distinct_lambdas = 385;
//		int number_of_distinct_lambdas = 100;
//		int number_of_distinct_lambdas = 400;
//		int number_of_distinct_lambdas = 100;
		int number_of_distinct_lambdas = 1000;

		int prefix_index = c_commoon_path.lastIndexOf("/");
		String prefix = c_commoon_path.substring(prefix_index + 1);

		String c_input_file_name = c_commoon_path + "/" + input_file_name;
		BufferedWriter bw = new BufferedWriter(new FileWriter(c_commoon_path + "/" + prefix + "__numLambdas_"
				+ (number_of_distinct_lambdas + 1) + "__" + input_file_name + "__" + all_c_input_file_name__set_T + "__"
				+ all_c_input_file_name__set_D + "__" + System.currentTimeMillis() + ".tsv"));

		double lambda_0 = 0.;
		double lambda_1 = 0.;
		double lambda_2 = 0.;
//		int T = 1000000000;
		int T = 10000;
//		int T = 1000;
//		int T = 100;
//		int T = 10;
//		int T = 1;
		double min_appx_factor = 0.999;
//		double min_appx_factor = 1.;

		String result = "";

		double[] all_lambda_1 = null;
		double[] all_lambda_2 = null;

		Set<String> all_solutions_found_so_far = new HashSet<String>();

		SuperGreedyPP algo = new SuperGreedyPP(c_input_file_name);
		algo.createGraph();
		algo.createSetsOfNodesToConsider(c_commoon_path + "/" + all_c_input_file_name__set_T,
				c_commoon_path + "/" + all_c_input_file_name__set_D, lambda_1, lambda_2);
		//
		double[] all_lambdas = SuperGreedyPP.get_max_ALL_lambdas(c_commoon_path, input_file_name,
				all_c_input_file_name__set_T, all_c_input_file_name__set_D, number_of_distinct_lambdas);
		//
		all_lambda_1 = all_lambdas;
		all_lambda_2 = all_lambdas;
		System.out.println("all_lambdas: " + Arrays.toString(all_lambdas));
		//
		lambda_0 = 1.;
		//
		for (int lamda_1_index = 0; lamda_1_index < all_lambda_1.length; lamda_1_index++) {
			//
			lambda_1 = all_lambda_1[lamda_1_index];
			//
			for (int lamda_2_index = 0; lamda_2_index < all_lambda_2.length; lamda_2_index++) {
				//
				lambda_2 = all_lambda_2[lamda_2_index];
				//
				if (lambda_0 != 1 && lambda_1 == 0 && lambda_2 == 0) {
					continue;
				}
				//
				algo.lambda_0 = lambda_0;
				algo.lambda_1 = lambda_1;
				algo.lambda_2 = lambda_2;
				//
				System.out.println(algo.lambda_0 + "\t" + algo.lambda_1 + "\t" + algo.lambda_2);
				//
				result = algo.runSuperGreedyPlusPlusOnlySolutionValue_2(T, min_appx_factor);
				int last_index = result.lastIndexOf("\t");
				String c_solution = result.substring(last_index + 1);
				if (all_solutions_found_so_far.add(c_solution)) {
					System.out.println("Numer of all_solutions_found_so_far: " + all_solutions_found_so_far.size());
					result = "";
					result += "DITH" + "\t";
					result += input_file_name + "\t";
					result += "dummy" + "\t";
					result += "dummy" + "\t";
					result += algo.lambda_1 + "\t";
					result += algo.lambda_2 + "\t";
					result += algo.map__field__value.get("AvgDeg(Sol)") + "\t";
					result += algo.map__field__value.get("AvgProx(Sol)") + "\t";
					result += algo.map__field__value.get("AvgDist(Sol)") + "\t";
					result += algo.map__field__value.get("OF(Sol)") + "\t";
					result += algo.map__field__value.get("UB") + "\t";
					result += algo.map__field__value.get("AppxFactor") + "\t";
					result += algo.map__field__value.get("t") + "\t";
					result += algo.map__field__value.get("ExecTimeMilliSec") + "\t";
					result += algo.map__field__value.get("T") + "\t";
					result += algo.map__field__value.get("gamma") + "\t";
					result += algo.map__field__value.get("|Sol|") + "\t";
					result += algo.map__field__value.get("Sol");
					bw.write(result);
					bw.write("\n");
					bw.flush();
				}
				//
			}
		}
		bw.close();
		return;
	}

	protected static void exp_SUPER_GREEDY_PP_ON_BLADE_with_lambdas_1_1_1(String c_commoon_path, String input_file_name,
			String all_c_input_file_name__set_T, String all_c_input_file_name__set_D) throws Exception {

		int prefix_index = c_commoon_path.lastIndexOf("/");
		String prefix = c_commoon_path.substring(prefix_index + 1);

		String c_input_file_name = c_commoon_path + "/" + input_file_name;
		BufferedWriter bw = new BufferedWriter(new FileWriter(
				c_commoon_path + "/" + prefix + "__" + input_file_name + "__" + all_c_input_file_name__set_T + "__"
						+ all_c_input_file_name__set_D + "__" + System.currentTimeMillis() + ".tsv"));

//		double lambda_0 = 0.;
		double lambda_1 = 0.;
		double lambda_2 = 0.;
//		int T = 1000000000;
		int T = 10000;
//		int T = 1000;
//		int T = 100;
//		int T = 10;
//		int T = 1;
		double min_appx_factor = 0.999;
//		double min_appx_factor = 1.;

		String result = "";

		Set<String> all_solutions_found_so_far = new HashSet<String>();

		SuperGreedyPP neg_dsg_c = new SuperGreedyPP(c_input_file_name);
		neg_dsg_c.createGraph();
		neg_dsg_c.createSetsOfNodesToConsider(c_commoon_path + "/" + all_c_input_file_name__set_T,
				c_commoon_path + "/" + all_c_input_file_name__set_D, lambda_1, lambda_2);
		// //
		//
		neg_dsg_c.lambda_0 = 1;
		neg_dsg_c.lambda_1 = 1;
		neg_dsg_c.lambda_2 = 1;
		//
		System.out.println("-_-");
		result = neg_dsg_c.runSuperGreedyPlusPlusOnlySolutionValue_2(T, min_appx_factor);
		System.out.println("_-_");
		int last_index = result.lastIndexOf("\t");
		String c_solution = result.substring(last_index + 1);
		if (all_solutions_found_so_far.add(c_solution)) {
			System.out.println("Numer of all_solutions_found_so_far: " + all_solutions_found_so_far.size());
			bw.write(result);
			bw.write("\n");
			bw.flush();
		}

		bw.close();
		return;
	}

	protected static void exp_RANDOM_SEARCH_ON_BLADE(String c_commoon_path, String input_file_name,
			String all_c_input_file_name__set_T, String all_c_input_file_name__set_D) throws Exception {

		int prefix_index = c_commoon_path.lastIndexOf("/");
		String prefix = c_commoon_path.substring(prefix_index + 1);

		String c_input_file_name = c_commoon_path + "/" + input_file_name;
		BufferedWriter bw = new BufferedWriter(new FileWriter(
				c_commoon_path + "/" + prefix + "__" + input_file_name + "__" + all_c_input_file_name__set_T + "__"
						+ all_c_input_file_name__set_D + "__" + System.currentTimeMillis() + ".tsv"));

		double lambda_0 = 0.;
		double lambda_1 = 0.;
		double lambda_2 = 0.;
//		int T = 1000000000;
		int T = 10000;
//		int T = 1000;
//		int T = 100;
//		int T = 10;
//		int T = 1;
		double min_appx_factor = 0.999;
//		double min_appx_factor = 1.;

		String result = "";
		//
		Set<String> all_solutions_found_so_far = new HashSet<String>();
		//
		SuperGreedyPP neg_dsg_c = new SuperGreedyPP(c_input_file_name);
		neg_dsg_c.createGraph();
		neg_dsg_c.createSetsOfNodesToConsider(c_commoon_path + "/" + all_c_input_file_name__set_T,
				c_commoon_path + "/" + all_c_input_file_name__set_D, lambda_1, lambda_2);
		//
		lambda_0 = 1.;
		//
		Random r = new Random();
		int lamda_1_MAX = 4539;
		int lamda_2_MAX = 3242;
		HashSet<Double> set_of_configurations = new HashSet<Double>(lamda_1_MAX * 2);
		while (set_of_configurations.size() < lamda_1_MAX * lamda_2_MAX) {
			//
			lambda_1 = r.nextInt(lamda_1_MAX + 1);
			lambda_2 = r.nextInt(lamda_2_MAX + 1);
			//
			if (!set_of_configurations.add(lambda_1 * lamda_2_MAX + lambda_2)) {
				continue;
			}
			//
			neg_dsg_c.lambda_0 = lambda_0;
			neg_dsg_c.lambda_1 = lambda_1;
			neg_dsg_c.lambda_2 = lambda_2;
			//
			System.out.println(neg_dsg_c.lambda_0 + "\t" + neg_dsg_c.lambda_1 + "\t" + neg_dsg_c.lambda_2);
			//
			result = neg_dsg_c.runSuperGreedyPlusPlusOnlySolutionValue_2(T, min_appx_factor);
			int last_index = result.lastIndexOf("\t");
			String c_solution = result.substring(last_index + 1);
			if (all_solutions_found_so_far.add(c_solution)) {
				System.out.println("Numer of all_solutions_found_so_far: " + all_solutions_found_so_far.size());
				System.out.println(result);
				bw.write(result);
				bw.write("\n");
				bw.flush();
			}
			//
		}
		bw.close();
		return;
	}

	protected int get_node_id_with_max_weight(double[] map__node_id__weight) {
		double c_max = -1.;
		int node_id_with_max_value = -1;
		double c_val = 0;
		for (int node_id = 0; node_id < map__node_id__weight.length; node_id++) {
			c_val = map__node_id__weight[node_id];
			if (c_val > c_max) {
				c_max = c_val;
				node_id_with_max_value = node_id;
			}
		}
		return node_id_with_max_value;
	}

	protected Set<Integer> get_ALL_nodes_id_with_max_weight(double[] map__node_id__weight) {
		Set<Integer> ALL_nodes_id_with_max_weight = new HashSet<Integer>();
		double c_max = -1.;
		double c_val = -1;
		for (int node_id = 0; node_id < map__node_id__weight.length; node_id++) {
			c_val = map__node_id__weight[node_id];
			if (c_val > c_max) {
				c_max = c_val;
			}
		}
		for (int node_id = 0; node_id < map__node_id__weight.length; node_id++) {
			c_val = map__node_id__weight[node_id];
			if (c_val == c_max) {
				ALL_nodes_id_with_max_weight.add(node_id);
			}
		}
		return ALL_nodes_id_with_max_weight;
	}

	protected boolean[] get_ego_network(int node_id) {
		Set<Integer> set_of_node_ids = new HashSet<Integer>(5);
		set_of_node_ids.add(node_id);
		return this.get_ego_network(set_of_node_ids);
	}

	protected boolean[] get_ego_network(Set<Integer> set_of_node_ids) {
		boolean[] subgraph_as_map__node_id__flag = new boolean[this.max_node_id + 1];
		Arrays.fill(subgraph_as_map__node_id__flag, false);
		//
		int[] node_id_adjacency_list = null;
		int c_neig_id = -1;
		int i = -1;
		for (int c_node_id : set_of_node_ids) {
			subgraph_as_map__node_id__flag[c_node_id] = true;
			//
			node_id_adjacency_list = this.graph_as_adjacency_list[c_node_id];
			c_neig_id = -1;
			for (i = 0; i < node_id_adjacency_list.length; i++) {
				c_neig_id = node_id_adjacency_list[i];
				//
				subgraph_as_map__node_id__flag[c_neig_id] = true;
				//
			}
		}
		//
		return subgraph_as_map__node_id__flag;
	}

	protected String run_on_a_SUBGRAPH_SuperGreedyPlusPlusOnlySolutionValue_2(boolean[] map__node_id__deleted_node)
			throws Exception {
		//
		long t_0 = System.currentTimeMillis();
		SuperGreedyPP sub_graph = this.createSubGraph(map__node_id__deleted_node);
		int[] map__sub_graph_node_id__node_id = new int[sub_graph.max_node_id + 1];
		int sub_graph_node_id;
		int graph_node_id;
		for (graph_node_id = 0; graph_node_id < this.map__old_node_id__new_node_id.length; graph_node_id++) {
			sub_graph_node_id = this.map__old_node_id__new_node_id[graph_node_id];
			if (sub_graph_node_id < 0) {
				continue;
			}
			map__sub_graph_node_id__node_id[sub_graph_node_id] = graph_node_id;
		}
		//
		sub_graph.lambda_0 = 1.;
		sub_graph.lambda_1 = 0.;
		sub_graph.lambda_2 = 0.;
		int T = 10000;
		double min_appx_factor = 0.99;
		sub_graph.runSuperGreedyPlusPlusOnlySolutionValue_2(T, min_appx_factor);
		//
		this.map__field__value = sub_graph.map__field__value;
		//
//		System.out.println(Arrays.toString(sub_graph.output_best_V));
		boolean[] c_solution = new boolean[this.max_node_id + 1];
		Arrays.fill(c_solution, false);
//		System.out.println("oiutyh sub_graph.max_node_id   " + sub_graph.max_node_id);
//		System.out.println("oiutyh sub_graph.E_w           " + sub_graph.E_w);
//		System.out.println("oiutyh sub_graph.output_best_V " + Arrays.toString(sub_graph.output_best_V));
		if (sub_graph.output_best_V != null) {
			for (int c_subgraph_node_id = 0; c_subgraph_node_id < sub_graph.output_best_V.length; c_subgraph_node_id++) {
				if (sub_graph.output_best_V[c_subgraph_node_id]) {
					graph_node_id = map__sub_graph_node_id__node_id[c_subgraph_node_id];
					c_solution[graph_node_id] = true;
				}
			}
		} else {
			Arrays.fill(c_solution, true);
		}
		long t_1 = System.currentTimeMillis();
		/// ofd
//		String result = "";
//		result = "";
//		result += "baseline" + "\t";
//		result += "input_file_name" + "\t";
//		result += "A" + "\t";
//		result += "B" + "\t";
//		result += sub_graph.lambda_1 + "\t";
//		result += sub_graph.lambda_2 + "\t";
//		result += sub_graph.map__field__value.get("AvgDeg(Sol)") + "\t";
//		result += sub_graph.map__field__value.get("AvgProx(Sol)") + "\t";
//		result += sub_graph.map__field__value.get("AvgDist(Sol)") + "\t";
//		result += sub_graph.map__field__value.get("OF(Sol)") + "\t";
//		result += sub_graph.map__field__value.get("UB") + "\t";
//		result += sub_graph.map__field__value.get("AppxFactor") + "\t";
//		result += sub_graph.map__field__value.get("t") + "\t";
//		result += sub_graph.map__field__value.get("ExecTimeMilliSec") + "\t";
//		result += sub_graph.map__field__value.get("T") + "\t";
//		result += sub_graph.map__field__value.get("gamma") + "\t";
//		result += sub_graph.map__field__value.get("|Sol|") + "\t";
//		result += sub_graph.map__field__value.get("Sol");
//		System.out.println("lkjhgfdsa");
//		System.out.println(result);
//		System.out.println("lkjhgfdsa");
		///

		int c_solution_size = computeSubgraphNumberOfNodes(c_solution);
		double density = computeSubgraphDensity(c_solution);
		double proximity_to_tau = computeSubgraphProximityToTau(c_solution);
		double distance_from_danger = computeSubgraphDistanceFromDanger(c_solution);
		double c_solution_value = computeSubgraphObjectiveFunctionValue(c_solution);
		String subgraph_as_string = convertSubgraphInSring3(c_solution);
		//
		this.map__field__value.put("AvgDeg(Sol)", "" + density);
		this.map__field__value.put("AvgProx(Sol)", "" + proximity_to_tau);
		this.map__field__value.put("AvgDist(Sol)", "" + distance_from_danger);
		this.map__field__value.put("OF(Sol)", "" + c_solution_value);
		this.map__field__value.put("UB", "-1");
		this.map__field__value.put("AppxFactor", "-1");
		// this.map__field__value.put("t", "" + t);
		this.map__field__value.put("ExecTimeMilliSec", "" + "" + (t_1 - t_0));
		// this.map__field__value.put("T", "" + T);
		this.map__field__value.put("gamma", "-1");
		this.map__field__value.put("|Sol|", "" + c_solution_size);
		this.map__field__value.put("Sol", subgraph_as_string);

		//
		String res = "";
		res = "-1" + "\t" + "-1" + "\t" + "-1" + "\t" + density + "\t" + proximity_to_tau + "\t" + distance_from_danger
				+ "\t" + c_solution_size + "\t" + "-1" + "\t" + "-1" + "\t" + "-1" + "\t" + "-1" + "\t" + (t_1 - t_0)
				+ "\t" + subgraph_as_string;
		//
		return res;
	}

	/**
	 * 
	 */
	protected String compute_baseline_MAX_DISTANCE_THEN_DENSITY() throws Exception {
		//
		int node_with_max_proximity = this.get_node_id_with_max_weight(this.map__node__weight_2);
		boolean[] map__node_id__deleted_node = this.get_ego_network(node_with_max_proximity);
		for (int i = 0; i < map__node_id__deleted_node.length; i++) {
			map__node_id__deleted_node[i] = !map__node_id__deleted_node[i];
		}
		//
//		System.out.println("node_with_max_proximity " + node_with_max_proximity);
//		System.out.println("node_with_max_proximity " + this.map__inner_id__outer_id[node_with_max_proximity]);
		//
		String res = this.run_on_a_SUBGRAPH_SuperGreedyPlusPlusOnlySolutionValue_2(map__node_id__deleted_node);
		//
		return res;
	}

	protected String compute_baseline_MAX_PROXIMITY_THEN_DENSITY() throws Exception {
		//
		Set<Integer> all_nodes_id_with_max_proximity = this.get_ALL_nodes_id_with_max_weight(this.map__node__weight_1);
		boolean[] map__node_id__deleted_node = this.get_ego_network(all_nodes_id_with_max_proximity);
		for (int i = 0; i < map__node_id__deleted_node.length; i++) {
			map__node_id__deleted_node[i] = !map__node_id__deleted_node[i];
		}
//		System.out.println(all_nodes_id_with_max_proximity);// ofd
//		for (int ssss : all_nodes_id_with_max_proximity) {
//			System.out.println(this.map__inner_id__outer_id[ssss]);
//		}
		//
//		System.out.println("node_with_max_distance " + node_with_max_distance);
//		System.out.println("node_with_max_distance " + this.map__inner_id__outer_id[node_with_max_distance]);
		//
		String res = this.run_on_a_SUBGRAPH_SuperGreedyPlusPlusOnlySolutionValue_2(map__node_id__deleted_node);
		//
		return res;
	}

	protected boolean are_equal(int node_id_i, int node_id_j) {
		boolean cond_1 = (this.map__node__weight_1[node_id_i] == this.map__node__weight_1[node_id_j]);
		boolean cond_2 = (this.map__node__weight_2[node_id_i] == this.map__node__weight_2[node_id_j]);
		return (cond_1 && cond_2);
	}

	protected boolean dominates(int node_id_i, int node_id_j) {
		if (are_equal(node_id_i, node_id_j)) {
			return false;
		}
		//
		double i_x = this.map__node__weight_1[node_id_i];
		double j_x = this.map__node__weight_1[node_id_j];
		double i_y = this.map__node__weight_2[node_id_i];
		double j_y = this.map__node__weight_2[node_id_j];
		//
		if ((i_x < j_x) || (i_y < j_y)) {
			return false;
		}
		//
		return true;
	}

	protected boolean[] skln() {
		boolean[] skln_as_map__node_id__flag = new boolean[this.max_node_id + 1];
		//
		Set<Integer> set__points_to_delete = new HashSet<Integer>();
		int node_id_i, node_id_j;
		for (node_id_i = 0; node_id_i <= this.max_node_id; node_id_i++) {
			for (node_id_j = node_id_i + 1; node_id_j <= this.max_node_id; node_id_j++) {
				//
				if (this.are_equal(node_id_i, node_id_j)) {
					continue;
				}
				if (this.dominates(node_id_i, node_id_j)) {
					set__points_to_delete.add(node_id_j);
				} else if (this.dominates(node_id_j, node_id_i)) {
					set__points_to_delete.add(node_id_i);
				}
				//
			}
		}
		//
//		int skln_size = 0;
		for (int node_id = 0; node_id < skln_as_map__node_id__flag.length; node_id++) {
			skln_as_map__node_id__flag[node_id] = true;
			if (set__points_to_delete.contains(node_id)) {
				skln_as_map__node_id__flag[node_id] = false;
			} else {
//				skln_size++;
//				System.out.println("skln elem: " + this.map__inner_id__outer_id[node_id] + "\t"
//						+ this.map__node__weight_1[node_id] + "\t" + this.map__node__weight_2[node_id]);
			}
		}
//		System.out.println("ewf skln_size=" + skln_size);
		//
		return skln_as_map__node_id__flag;
	}

	protected String compute_baseline_MAX_PROXIMITY_AND_DISTANCE_THEN_DENSITY() throws Exception {
		// skln computation.
		boolean[] map__node_id__deleted_node = this.skln();
		for (int i = 0; i < map__node_id__deleted_node.length; i++) {
			map__node_id__deleted_node[i] = !map__node_id__deleted_node[i];
		}
		//
		String res = this.run_on_a_SUBGRAPH_SuperGreedyPlusPlusOnlySolutionValue_2(map__node_id__deleted_node);
		//
		return res;
	}

	protected String compute_baseline_MAX_DENSITY() throws Exception {
		//
		double original_lambda_0 = this.lambda_0;
		double original_lambda_1 = this.lambda_1;
		double original_lambda_2 = this.lambda_2;
		//
		this.lambda_0 = 1.;
		this.lambda_1 = 0.;
		this.lambda_2 = 0.;
		//
		int T = 10000;
		double min_appx_factor = 0.99;
		//
		String result = this.runSuperGreedyPlusPlusOnlySolutionValue_2(T, min_appx_factor);
		//
		this.lambda_0 = original_lambda_0;
		this.lambda_1 = original_lambda_1;
		this.lambda_2 = original_lambda_2;
		//
		double c_solution_value = computeSubgraphObjectiveFunctionValue(this.output_best_V);
		//
//		this.map__field__value.put("AvgDeg(Sol)", "" + density);
//		this.map__field__value.put("AvgProx(Sol)", "" + proximity_to_tau);
//		this.map__field__value.put("AvgDist(Sol)", "" + distance_from_danger);
		this.map__field__value.put("OF(Sol)", "" + c_solution_value);
//		this.map__field__value.put("UB", "-1");
//		this.map__field__value.put("AppxFactor", "-1");
		// this.map__field__value.put("t", "" + t);
//		this.map__field__value.put("ExecTimeMilliSec", "" + "" + (t_1 - t_0));
		// this.map__field__value.put("T", "" + T);
//		this.map__field__value.put("gamma", "-1");
//		this.map__field__value.put("|Sol|", "" + c_solution_size);
//		this.map__field__value.put("Sol", subgraph_as_string);

		//
		return result;
	}

	protected String compute_baseline_MAX_DISTANCE() throws Exception {
		long t_0 = System.currentTimeMillis();
		int node_with_max_proximity = this.get_node_id_with_max_weight(this.map__node__weight_1);
		boolean[] c_solution = new boolean[this.max_node_id + 1];
		Arrays.fill(c_solution, false);
		c_solution[node_with_max_proximity] = true;
		long t_1 = System.currentTimeMillis();
		//
//		System.out.println("node_with_max_proximity " + node_with_max_proximity);
//		System.out.println("node_with_max_proximity " + this.map__inner_id__outer_id[node_with_max_proximity]);
		//
		double c_solution_size = computeSubgraphNumberOfNodes(c_solution);
		double density = computeSubgraphDensity(c_solution);
		double proximity_to_tau = computeSubgraphProximityToTau(c_solution);
		double distance_from_danger = computeSubgraphDistanceFromDanger(c_solution);
		String subgraph_as_string = convertSubgraphInSring3(c_solution);
		//
		String res = "";
		res = "-1" + "\t" + "-1" + "\t" + "-1" + "\t" + density + "\t" + proximity_to_tau + "\t" + distance_from_danger
				+ "\t" + c_solution_size + "\t" + "-1" + "\t" + "-1" + "\t" + "-1" + "\t" + "-1" + "\t" + (t_1 - t_0)
				+ "\t" + subgraph_as_string;
		//
		return res;
	}

	protected String compute_baseline_MAX_PROXIMITY() throws Exception {
		long t_0 = System.currentTimeMillis();
		int node_with_max_distance = this.get_node_id_with_max_weight(this.map__node__weight_2);
		boolean[] c_solution = new boolean[this.max_node_id + 1];
		Arrays.fill(c_solution, false);
		c_solution[node_with_max_distance] = true;
		long t_1 = System.currentTimeMillis();
		//
//		System.out.println("node_with_max_proximity " + node_with_max_distance);
//		System.out.println("node_with_max_proximity " + this.map__inner_id__outer_id[node_with_max_distance]);
		//
		double c_solution_size = computeSubgraphNumberOfNodes(c_solution);
		double density = computeSubgraphDensity(c_solution);
		double proximity_to_tau = computeSubgraphProximityToTau(c_solution);
		double distance_from_danger = computeSubgraphDistanceFromDanger(c_solution);
		String subgraph_as_string = convertSubgraphInSring3(c_solution);
		//
		String res = "";
		res = "-1" + "\t" + "-1" + "\t" + "-1" + "\t" + density + "\t" + proximity_to_tau + "\t" + distance_from_danger
				+ "\t" + c_solution_size + "\t" + "-1" + "\t" + "-1" + "\t" + "-1" + "\t" + "-1" + "\t" + (t_1 - t_0)
				+ "\t" + subgraph_as_string;
		//
		return res;
	}

	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////

	protected static void exp_GRID_SEARCH(String algo) throws Exception {
		String c_commoon_path = "";
		String all_c_input_file_name__set_T = "";
		String all_c_input_file_name__set_D = "";
		String input_file_name = "";

		/////////////////////////////////////////////////////////////////////////////////////////////
//		c_commoon_path = "/Users/ikki/Downloads/ISI/retweet_russia_march";
//		input_file_name = "russia_march.tsv";
//		all_c_input_file_name__set_T = c_commoon_path + "/" + "T__TSPR.tsv";
//		all_c_input_file_name__set_D = c_commoon_path + "/" + "D__1_minus_TSPR.tsv";
		/////////////////////////////////////////////////////////////////////////////////////////////

		/////////////////////////////////////////////////////////////////////////////////////////////
//		c_commoon_path = "/Users/ikki/Downloads/ISI/vaxnovax";
		c_commoon_path = "/home/fazzone/sofarsoclose/vaxnovax";
		input_file_name = "follow_graph.tsv";
////		all_c_input_file_name__set_T = c_commoon_path + "/" + "NOVAX_GT__MAX_minus_MSSP.tsv";
////		all_c_input_file_name__set_D = c_commoon_path + "/" + "PROVAX_GT__MSSP.tsv";
		all_c_input_file_name__set_T = c_commoon_path + "/" + "NOVAX_GT__TSPR.tsv";
		all_c_input_file_name__set_D = c_commoon_path + "/" + "PROVAX_GT__1_minus_TSPR.tsv";
		/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
//		c_commoon_path = "/Users/ikki/Downloads/ISI/vaxnovax/retweet";
//		c_commoon_path = "/home/fazzone/sofarsoclose/vaxnovax/retweet";
//		input_file_name = "retweet_edges.tsv";
////		all_c_input_file_name__set_T = c_commoon_path + "/" + "NOVAX_GT__MAX_minus_MSSP.tsv";
////		all_c_input_file_name__set_D = c_commoon_path + "/" + "PROVAX_GT__MSSP.tsv";
//		all_c_input_file_name__set_T = c_commoon_path + "/" + "NOVAX_GT__TSPR.tsv";
//		all_c_input_file_name__set_D = c_commoon_path + "/" + "PROVAX_GT__1_minus_TSPR.tsv";
/////////////////////////////////////////////////////////////////////////////////////////////

// DSG as value
// MAX-proximity as value
// MAX-distance as value
//		 * 
// DSG as solution
// DSG on the Tau ego-network as a solution.
// DSG on the ego-network of the the furthest node.
//		 * 
// Extraction of the best nodes in terms of both proximity and distance, then DSG on their induced subgraph. 

		String c_input_file_name = c_commoon_path + "/" + input_file_name;
		BufferedWriter bw = new BufferedWriter(
				new FileWriter(c_commoon_path + "/" + input_file_name + "__" + System.currentTimeMillis() + ".tsv"));

		double lambda_0 = 0.;
		double lambda_1 = 0.;
		double lambda_2 = 0.;
//		int T = 1000000000;
		int T = 10000;
//		int T = 1000;
//		int T = 100;
//		int T = 10;
//		int T = 1;
		double min_appx_factor = 0.999;
//		double min_appx_factor = 1.;

		String result = "";
//		double[] all_lambda_0 = { 1, 2, 3, 4, 5, 6, 7 };
//		double[] all_lambda_1 = { 0, 1, 2, 3, 4, 5, 6, 7 };
//		double[] all_lambda_2 = { 0, 1, 2, 3, 4, 5, 6, 7 };
//		double[] all_lambda_1 = { 0, 1, 10, 100, 1000, 10000, 100000 };
//		double[] all_lambda_2 = { 0, 1, 10, 100, 1000, 10000, 100000 };
//		double[] all_lambda_1 = { 0, 10000, 100000 };
//		double[] all_lambda_2 = { 0, 10000, 100000 };

//		double[] all_lambda_0 = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30};
		double[] all_lambda_0 = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400,
				500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000 };
		double[] all_lambda_1 = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400,
				500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000 };
		double[] all_lambda_2 = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400,
				500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000 };

		Set<String> all_solutions_found_so_far = new HashSet<String>();

		SuperGreedyPP neg_dsg_c = new SuperGreedyPP(c_input_file_name);
		neg_dsg_c.createGraph();
		neg_dsg_c.createSetsOfNodesToConsider(all_c_input_file_name__set_T, all_c_input_file_name__set_D, lambda_1,
				lambda_2);
		for (int lamda_0_index = 0; lamda_0_index < all_lambda_0.length; lamda_0_index++) {
//		for (lambda_0 = 0; lambda_0 <= 30; lambda_0++) {
			//
			lambda_0 = all_lambda_0[lamda_0_index];
			//
			for (int lamda_1_index = 0; lamda_1_index < all_lambda_1.length; lamda_1_index++) {
				//
				lambda_1 = all_lambda_1[lamda_1_index];
				//
				for (int lamda_2_index = 0; lamda_2_index < all_lambda_2.length; lamda_2_index++) {
					//
					lambda_2 = all_lambda_2[lamda_2_index];
					//
					if (lambda_0 != 1 && lambda_1 == 0 && lambda_2 == 0) {
						continue;
					}
					//
					neg_dsg_c.lambda_0 = lambda_0;
					neg_dsg_c.lambda_1 = lambda_1;
					neg_dsg_c.lambda_2 = lambda_2;
					//
					System.out.println(neg_dsg_c.lambda_0 + "\t" + neg_dsg_c.lambda_1 + "\t" + neg_dsg_c.lambda_2);
					//
					result = neg_dsg_c.runSuperGreedyPlusPlusOnlySolutionValue_2(T, min_appx_factor);
					int last_index = result.lastIndexOf("\t");
					String c_solution = result.substring(last_index + 1);
					if (all_solutions_found_so_far.add(c_solution)) {
						System.out.println("Numer of all_solutions_found_so_far: " + all_solutions_found_so_far.size());
						bw.write(result);
						bw.write("\n");
						bw.flush();
					}
					//
				}
			}
		}
		bw.close();
		return;
	}

	protected static void exp_GRID_SEARCH_ON_BLADE(String c_commoon_path, String input_file_name,
			String all_c_input_file_name__set_T, String all_c_input_file_name__set_D) throws Exception {

		int prefix_index = c_commoon_path.lastIndexOf("/");
		String prefix = c_commoon_path.substring(prefix_index + 1);

		String c_input_file_name = c_commoon_path + "/" + input_file_name;
		BufferedWriter bw = new BufferedWriter(new FileWriter(
				c_commoon_path + "/" + prefix + "__" + input_file_name + "__" + all_c_input_file_name__set_T + "__"
						+ all_c_input_file_name__set_D + "__" + System.currentTimeMillis() + ".tsv"));

		double lambda_0 = 0.;
		double lambda_1 = 0.;
		double lambda_2 = 0.;
//		int T = 1000000000;
		int T = 10000;
//		int T = 1000;
//		int T = 100;
//		int T = 10;
//		int T = 1;
		double min_appx_factor = 0.999;
//		double min_appx_factor = 1.;

		String result = "";
//		double[] all_lambda_0 = { 1, 2, 3, 4, 5, 6, 7 };
//		double[] all_lambda_1 = { 0, 1, 2, 3, 4, 5, 6, 7 };
//		double[] all_lambda_2 = { 0, 1, 2, 3, 4, 5, 6, 7 };
//		double[] all_lambda_1 = { 0, 1, 10, 100, 1000, 10000, 100000 };
//		double[] all_lambda_2 = { 0, 1, 10, 100, 1000, 10000, 100000 };
//		double[] all_lambda_1 = { 0, 10000, 100000 };
//		double[] all_lambda_2 = { 0, 10000, 100000 };

//		double[] all_lambda_0 = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30};
//		double[] all_lambda_0 = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400,
//				500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000 };
//		double[] all_lambda_1 = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400,
//				500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000 };
//		double[] all_lambda_2 = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400,
//				500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000 };

//		double[] all_lambda_0 = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2, 3, 4, 5 };
//		double[] all_lambda_1 = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };
//		double[] all_lambda_2 = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100 };

//		double[] all_lambda_0 = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400,
//				500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000 };
//		double[] all_lambda_1 = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400,
//				500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000 };
//		double[] all_lambda_2 = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400,
//				500, 600, 700, 800, 900, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000 };

		double[] all_lambda_0 = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400,
				500, 600, 700, 800, 900, 1000 };
		double[] all_lambda_1 = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400,
				500, 600, 700, 800, 900, 1000 };
		double[] all_lambda_2 = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400,
				500, 600, 700, 800, 900, 1000 };

		Set<String> all_solutions_found_so_far = new HashSet<String>();

		SuperGreedyPP neg_dsg_c = new SuperGreedyPP(c_input_file_name);
		neg_dsg_c.createGraph();
		neg_dsg_c.createSetsOfNodesToConsider(c_commoon_path + "/" + all_c_input_file_name__set_T,
				c_commoon_path + "/" + all_c_input_file_name__set_D, lambda_1, lambda_2);
		for (int lamda_0_index = 0; lamda_0_index < all_lambda_0.length; lamda_0_index++) {
			//
			lambda_0 = all_lambda_0[lamda_0_index];
			//
			for (int lamda_1_index = 0; lamda_1_index < all_lambda_1.length; lamda_1_index++) {
				//
				lambda_1 = all_lambda_1[lamda_1_index];
				//
				for (int lamda_2_index = 0; lamda_2_index < all_lambda_2.length; lamda_2_index++) {
					//
					lambda_2 = all_lambda_2[lamda_2_index];
					//
					if (lambda_0 != 1 && lambda_1 == 0 && lambda_2 == 0) {
						continue;
					}
					//
					neg_dsg_c.lambda_0 = lambda_0;
					neg_dsg_c.lambda_1 = lambda_1;
					neg_dsg_c.lambda_2 = lambda_2;
					//
					System.out.println(neg_dsg_c.lambda_0 + "\t" + neg_dsg_c.lambda_1 + "\t" + neg_dsg_c.lambda_2);
					//
					result = neg_dsg_c.runSuperGreedyPlusPlusOnlySolutionValue_2(T, min_appx_factor);
					int last_index = result.lastIndexOf("\t");
					String c_solution = result.substring(last_index + 1);
					if (all_solutions_found_so_far.add(c_solution)) {
						System.out.println("Numer of all_solutions_found_so_far: " + all_solutions_found_so_far.size());
						bw.write(result);
						bw.write("\n");
						bw.flush();
					}
					//
				}
			}
		}
		bw.close();
		return;
	}

	public String runSuperGreedyPlusPlusOnlySolutionValue_2(int T, double min_appx_factor) throws Exception {
		return this.runSuperGreedyPlusPlusOnlySolutionValue_2(T, min_appx_factor, false);
	}

	public String runSuperGreedyPlusPlusOnlySolutionValue_2(int T, double min_appx_factor, boolean verbose)
			throws Exception {
		long t_0 = 0;
		long t_1 = 0;
		long t_x = 0;
		t_0 = System.currentTimeMillis();
		String out_RESULT_string = "";
		String out_string = "";
		// keep
		double[] map__node_id__load = new double[this.max_node_id + 1];
		//
		//
		double solution_quality = 0.;
		double c_node_load;
		double c_solution_value = 0.;
		int c_solution_size = 0;
		double best_solution_so_far_value = 0.;
		int best_solution_so_far_size = 0;
		double UB = Double.MAX_VALUE;
		double c_iteration_MAX_LOAD = Double.MIN_VALUE;
		double c_iteration_UB = Double.MAX_VALUE;
		boolean last_iteration = false;
		boolean just_after_the_last_iteration = false;
		int t;
		for (t = 1; t <= T; t++) {
			if (just_after_the_last_iteration) {
				break;
			}
			//
			this.singleIterationOfSuperGreedyPlusPlusOnlySolutionValue(map__node_id__load);
			c_solution_value = this.output_best_solution_value;
			c_solution_size = this.output_best_V_size;
			//
			// Update best solution encountered so far.
			if (c_solution_value > best_solution_so_far_value) {
				best_solution_so_far_value = c_solution_value;
				best_solution_so_far_size = c_solution_size;
				this.computeSolutionRepresentation(this.index_for_best_solution_set_of_nodes, this.list__removed_nodes);
			}
			//
			// compute LB and UB
			c_iteration_MAX_LOAD = -1;
			for (int c_node_id = 0; c_node_id < map__node_id__load.length; c_node_id++) {
				//
				c_node_load = map__node_id__load[c_node_id];
				//
				c_iteration_MAX_LOAD = (c_node_load > c_iteration_MAX_LOAD ? c_node_load : c_iteration_MAX_LOAD);
				//
			}
			//
			c_iteration_UB = c_iteration_MAX_LOAD / t;
			UB = (c_iteration_UB < UB ? c_iteration_UB : UB);
			//
			solution_quality = best_solution_so_far_value / UB;
			//
			last_iteration = (t == T) || (solution_quality >= min_appx_factor);
			//
			if (last_iteration) {
				if (this.output_best_V == null) {
//					System.out.println("rjkwnfnksdl " + this.output_best_V + "\t" + this.lambda_0 + "\t" + this.lambda_0
//							+ "\t" + this.lambda_0 + "\t");
					this.output_best_V = new boolean[this.max_node_id + 1];
					Arrays.fill(this.output_best_V, true);
				}
				t_1 = System.currentTimeMillis();
				just_after_the_last_iteration = true;
				//
				double density = computeSubgraphDensity(this.output_best_V);
				double proximity_to_tau = computeSubgraphProximityToTau(this.output_best_V);
				double distance_from_danger = computeSubgraphDistanceFromDanger(this.output_best_V);
				String subgraph_as_string = convertSubgraphInSring3(this.output_best_V);
				out_RESULT_string = this.lambda_0 + "\t" + this.lambda_1 + "\t" + this.lambda_2 + "\t" + density + "\t"
						+ proximity_to_tau + "\t" + distance_from_danger + "\t" + best_solution_so_far_size + "\t"
						+ best_solution_so_far_value + "\t" + UB + "\t" + solution_quality + "\t" + t + "\t"
						+ (t_1 - t_0) + "\t" + subgraph_as_string;

				this.map__field__value = new HashMap<String, String>(30);
				this.map__field__value.put("AvgDeg(Sol)", "" + density);
				this.map__field__value.put("AvgProx(Sol)", "" + proximity_to_tau);
				this.map__field__value.put("AvgDist(Sol)", "" + distance_from_danger);
				this.map__field__value.put("OF(Sol)", "" + best_solution_so_far_value);
				this.map__field__value.put("UB", "" + UB);
				this.map__field__value.put("AppxFactor", "" + solution_quality);
				this.map__field__value.put("t", "" + t);
				this.map__field__value.put("ExecTimeMilliSec", "" + "" + (t_1 - t_0));
				this.map__field__value.put("T", "" + T);
				this.map__field__value.put("gamma", "" + min_appx_factor);
				this.map__field__value.put("|Sol|", "" + best_solution_so_far_size);
				this.map__field__value.put("Sol", subgraph_as_string);

//				
//				input:  [method, G, A, B, lambda_1, lambda_2, file_prox, file_dist]
//				output: [method, G, A, B, lambda_1, lambda_2, AvgDeg(Sol), AvgProx(Sol), AvgDist(Sol), OF(Sol), UB, AppxFactor, t, ExecTimeMilliSec, T, gamma, |Sol|, Sol]

//				System.out.println();
//				System.out.println(out_string);
			}
			//
			t_x = System.currentTimeMillis();
			//
			if (verbose) {
				if (last_iteration) {
					out_string = t + "\t" + (t_1 - t_0) + "\t" + best_solution_so_far_value + "\t" + UB + "\t"
							+ solution_quality + "\t" + best_solution_so_far_size;
				} else {
					out_string = t + "\t" + (t_x - t_0) + "\t" + best_solution_so_far_value + "\t" + UB + "\t"
							+ solution_quality + "\t" + best_solution_so_far_size;
				}
				System.out.println(out_string);
			}
			//
		}
		//
		return out_RESULT_string;
	}

	public boolean[] runSuperGreedyPlusPlusOnlySolutionValue(int T, double min_appx_factor) throws Exception {
		return runSuperGreedyPlusPlusOnlySolutionValue(T, min_appx_factor, null, -1);
	}

	public boolean[] runSuperGreedyPlusPlusOnlySolutionValue(int T, double min_appx_factor, boolean[] map__node__tau,
			int num_nodes_in_tau) throws Exception {
		// keep
		double[] map__node_id__load = new double[this.max_node_id + 1];
		//
		long t_0, t_1;
		//
		String out_string = "";
		out_string = "t" + "\t" + "best_solution_value_so_far" + "\t" + "UB" + "\t"
				+ "SIZE_of_best_solution_so_far(num_nodes)" + "\t"
				+ "Solution_Quality(BestSolValue/UP + \"\\t\" + \"exec_time msec\" )";
		System.out.println(out_string);
		//
		double solution_quality = 0.;
		double c_node_load;
		double c_solution_value = 0.;
		int c_solution_size = 0;
		double best_solution_so_far_value = 0.;
		int best_solution_so_far_size = 0;
		double UB = Double.MAX_VALUE;
		double c_iteration_MAX_LOAD = Double.MIN_VALUE;
		double c_iteration_UB = Double.MAX_VALUE;
		boolean last_iteration = false;
		int t;
		for (t = 1; t <= T; t++) {
			last_iteration = (t == T) || (solution_quality >= min_appx_factor);
			//
			t_0 = System.currentTimeMillis();
			if (map__node__tau == null) {
				this.singleIterationOfSuperGreedyPlusPlusOnlySolutionValue(map__node_id__load);
			} else {
				this.singleIterationOfSuperGreedyPlusPlusOnlySolutionValue(map__node_id__load, map__node__tau,
						num_nodes_in_tau);
			}
			t_1 = System.currentTimeMillis();
			c_solution_value = this.output_best_solution_value;
			c_solution_size = this.output_best_V_size;
			//
			// Update best solution encountered so far.
			if (c_solution_value > best_solution_so_far_value) {
				best_solution_so_far_value = c_solution_value;
				best_solution_so_far_size = c_solution_size;
				this.computeSolutionRepresentation(this.index_for_best_solution_set_of_nodes, this.list__removed_nodes);
			}
			//
			// compute LB and UB
			c_iteration_MAX_LOAD = -1;
			for (int c_node_id = 0; c_node_id < map__node_id__load.length; c_node_id++) {
				//
				c_node_load = map__node_id__load[c_node_id];
				//
				c_iteration_MAX_LOAD = (c_node_load > c_iteration_MAX_LOAD ? c_node_load : c_iteration_MAX_LOAD);
				//
			}
			//
			c_iteration_UB = c_iteration_MAX_LOAD / t;
			UB = (c_iteration_UB < UB ? c_iteration_UB : UB);
			//
			solution_quality = best_solution_so_far_value / UB;
			out_string = t + "\t" + best_solution_so_far_value + "\t" + UB + "\t" + best_solution_so_far_size + "\t"
					+ solution_quality + "\t" + (t_1 - t_0);
			System.out.println(out_string);
			if (last_iteration) {
				System.out.println(out_string);
//				System.out.println(Arrays.toString(this.output_best_V));
//				String complementingNodes = printComplementingNodes(this.output_best_V);
//				System.out.println();
//				System.out.println(complementingNodes);
//				BufferedWriter bw = new BufferedWriter(new FileWriter("/Users/ikki/Downloads/ISI/"
//						+ "retweet_graph_mothersday_threshold_largest_CC__complementingNodes__"
//						+ System.currentTimeMillis() + ".tsv"));
//				bw.write(complementingNodes);
//				bw.close();
				double[] map__class__num_nodes = this.computeSubgraphBalance(this.output_best_V);
				double density = computeSubgraphDensity(this.output_best_V);
				System.out.println("density = " + density);
				System.out.println("size2   = " + computeSubgraphNumNodes(this.output_best_V));
				System.out.println("lambda_1= " + this.lambda_1);
				System.out.println("lambda_2= " + this.lambda_2);
				System.out.println(Arrays.toString(map__class__num_nodes));
				String subgraph_as_string = this.convertSubgraphInSring(this.output_best_V);
				System.out.println(subgraph_as_string);
//				subgraph_as_string = this.convertSubgraphInSring2(this.output_best_V);
				subgraph_as_string = convertSubgraphInSring3(this.output_best_V);
				System.out.println(subgraph_as_string);
				System.out.println();
//				String aaa = "" + "A_num_nodes = 25  # [0, 25) dsg that belongs to the NORTH.\n"
//						+ "B_num_nodes = 25  # [25, 50)other community of the NORTH.\n"
//						+ "C_num_nodes = 100  # [50, 150) superset of Dengerous nodes, of course in the NORTH.\n" + "\n"
//						+ "D_num_nodes = 100  # [150, 250) superset of Ground-Truth nodes, of course in the SOUTH.\n"
//						+ "E_num_nodes = 25  # [250, 275) other community of the SOUTH, but not the best solution in terms of density.\n"
//						+ "F_num_nodes = 25  # [275, 299) BEST SOLUTION";
//				System.out.println(aaa);
				System.out.println();
				break;
			}
			//
		}
		//
		return null;
	}

	protected static void exp_ON_TWEETTER(String algo) throws Exception {
		// keep
		String c_commoon_path = "/Users/ikki/Downloads/ISI/retweet_russia_march";
		String all_c_input_file_name__set_T = "";
		String all_c_input_file_name__set_D = "";
		String[] all_c_input_file_name = { "russia_march.tsv" };
		//
		//
//		all_c_input_file_name__set_T = c_commoon_path + "/" + "T.tsv";
		all_c_input_file_name__set_T = c_commoon_path + "/" + "T__TSPR.tsv";
		all_c_input_file_name__set_D = c_commoon_path + "/" + "D__1_minus_TSPR.tsv";
//		all_c_input_file_name__set_D = c_commoon_path + "/" + "D__MSSP.tsv";
//		
//		all_c_input_file_name__set_T = c_commoon_path + "/" + "T__TSPR.tsv";
//		all_c_input_file_name__set_D = c_commoon_path + "/" + "D__1_minus_TSPR.tsv";
//		
//		all_c_input_file_name__set_T = c_commoon_path + "/" + "T__MAX_minus_MSSP.tsv";
//		all_c_input_file_name__set_D = c_commoon_path + "/" + "D__MSSP.tsv";
//
//		all_c_input_file_name__set_T = c_commoon_path + "/" + "T__MAX_minus_MSSP.tsv";
//		all_c_input_file_name__set_D = c_commoon_path + "/" + "D__1_minus_TSPR.tsv";
//
//		all_c_input_file_name__set_T = c_commoon_path + "/" + "T__MAX_minus_MSSP.tsv";
//		all_c_input_file_name__set_D = c_commoon_path + "/" + "D__MSSP.tsv";

		double lambda_1 = 10.;
		double lambda_2 = 10000.;

		/*
		 * double lambda_1 = 100.; double lambda_2 = 5000.;
		 * 
		 * 
		 * double lambda_1 = 100.; double lambda_2 = 100000 / 4.;
		 * 
		 * double lambda_1 = 0.; double lambda_2 = 100000 / 4.;
		 * 
		 * 
		 */
		String c_input_file_name;
		for (String input_file_name : all_c_input_file_name) {
			c_input_file_name = c_commoon_path + "/" + input_file_name;
			//
			//
			//
			System.out.println();
			System.out.println("-----------------------------------------");
			System.out.println(c_input_file_name);
			System.out.println();
			//
//			System.gc();
			SuperGreedyPP neg_dsg_c = new SuperGreedyPP(c_input_file_name);
			neg_dsg_c.createGraph();
			neg_dsg_c.createSetsOfNodesToConsider(all_c_input_file_name__set_T, all_c_input_file_name__set_D, lambda_1,
					lambda_2);
//			NegDSGCharikarBabis neg_dsg_c = loadSerializeGraphFromFile(c_input_file_name);
//			System.gc();
			//
			//
			//
			System.out.println("input_file_name: " + input_file_name);
			//

//			int T = 1000000000;
			int T = 10000;
//			int T = 1000;
//			int T = 100;
//			int T = 10;
			double min_appx_factor = 0.999;
//			double min_appx_factor = 1.;
			neg_dsg_c.runSuperGreedyPlusPlusOnlySolutionValue(T, min_appx_factor);
		}
		return;
	}

	protected static void exp_ON_MODEL(String algo) throws Exception {
		// keep
		String c_commoon_path = "/Users/ikki/Dropbox/NegDSG/sw/MODEL_3";
		String all_c_input_file_name__set_T = "";
		String all_c_input_file_name__set_D = "";
		String[] all_c_input_file_name = { "sbm.tsv" };
		//
		//
		all_c_input_file_name__set_T = c_commoon_path + "/" + "T.tsv";
		all_c_input_file_name__set_D = c_commoon_path + "/" + "D__MSSP.tsv";
//		
//		all_c_input_file_name__set_T = c_commoon_path + "/" + "T__TSPR.tsv";
//		all_c_input_file_name__set_D = c_commoon_path + "/" + "D__1_minus_TSPR.tsv";
//		
//		all_c_input_file_name__set_T = c_commoon_path + "/" + "T__MAX_minus_MSSP.tsv";
//		all_c_input_file_name__set_D = c_commoon_path + "/" + "D__MSSP.tsv";
//
//		all_c_input_file_name__set_T = c_commoon_path + "/" + "T__MAX_minus_MSSP.tsv";
//		all_c_input_file_name__set_D = c_commoon_path + "/" + "D__1_minus_TSPR.tsv";
//
//		all_c_input_file_name__set_T = c_commoon_path + "/" + "T__MAX_minus_MSSP.tsv";
//		all_c_input_file_name__set_D = c_commoon_path + "/" + "D__MSSP.tsv";

		double lambda_1 = 0.;
		double lambda_2 = 1.;

		String c_input_file_name;
		for (String input_file_name : all_c_input_file_name) {
			c_input_file_name = c_commoon_path + "/" + input_file_name;
			//
			//
			//
			System.out.println();
			System.out.println("-----------------------------------------");
			System.out.println(c_input_file_name);
			System.out.println();
			//
//			System.gc();
			SuperGreedyPP neg_dsg_c = new SuperGreedyPP(c_input_file_name);
			neg_dsg_c.createGraph();
			neg_dsg_c.createSetsOfNodesToConsider(all_c_input_file_name__set_T, all_c_input_file_name__set_D, lambda_1,
					lambda_2);
//			NegDSGCharikarBabis neg_dsg_c = loadSerializeGraphFromFile(c_input_file_name);
//			System.gc();
			//
			//
			//
			System.out.println("input_file_name: " + input_file_name);
			//

//			int T = 1000000000;
			int T = 10000;
//			int T = 1000;
//			int T = 100;
//			int T = 10;
			double min_appx_factor = 0.999;
//			double min_appx_factor = 1.;
			neg_dsg_c.runSuperGreedyPlusPlusOnlySolutionValue(T, min_appx_factor);
		}
		return;
	}

	protected static void exp_12(String algo) throws Exception {
		// keep
//		String c_commoon_path = "/Users/ikki/DATASETS/NegDSG/vaxnovax";
//		String[] all_c_input_file_name = { "follow_graph.tsv" };
//		String[] all_c_input_file_name = { "random_graph__10.tsv" };
		//
//		String all_c_input_file_name_NODES_TO_CONSIDER_SET_1 = c_commoon_path + "/" + "dummy.tsv";
//		String all_c_input_file_name_NODES_TO_CONSIDER_SET_1 = c_commoon_path + "/" + "random_graph__10__SET_T.tsv";
//		String all_c_input_file_name_NODES_TO_CONSIDER_SET_1 = c_commoon_path + "/" + "node_id__PROVAX_SCORE.tsv";
//		String all_c_input_file_name_NODES_TO_CONSIDER_SET_1 = c_commoon_path + "/" + "node_id__SELECTED_PROVAX_MAX_SCORE.tsv";
//		String all_c_input_file_name_NODES_TO_CONSIDER_SET_1 = c_commoon_path + "/" + "node_id__SELECTED_NOVAX_MAX_SCORE__2.tsv";
//		String all_c_input_file_name_NODES_TO_CONSIDER_SET_1 = c_commoon_path + "/" + "node_id__SELECTED_PROVAX_DISTANCE.tsv";
//		String all_c_input_file_name_NODES_TO_CONSIDER_SET_1 = c_commoon_path + "/" + "node_id__SELECTED_PROVAX_DISTANCE__QUAD.tsv";

//		String all_c_input_file_name_NODES_TO_CONSIDER_SET_2 = c_commoon_path + "/" + "dummy.tsv";
//		String all_c_input_file_name_NODES_TO_CONSIDER_SET_2 = c_commoon_path + "/"
//				+ "random_graph__10__DISTANCES_FROM_SET_D.tsv";
//		String all_c_input_file_name_NODES_TO_CONSIDER_SET_2 = c_commoon_path + "/" + "node_id__NOVAX_SCORE.tsv";
//		String all_c_input_file_name_NODES_TO_CONSIDER_SET_2 = c_commoon_path + "/" + "node_id__NOVAX_DISTANCE.tsv";
//		String all_c_input_file_name_NODES_TO_CONSIDER_SET_2 = c_commoon_path + "/" + "node_id__SELECTED_NOVAX_DISTANCE.tsv";
//		String all_c_input_file_name_NODES_TO_CONSIDER_SET_2 = c_commoon_path + "/" + "node_id__SELECTED_NOVAX_DISTANCE__QUAD.tsv";

		String c_commoon_path = "/Users/ikki/Downloads/ISI/retweet_russia_march/retweet";
		String all_c_input_file_name__set_T = "";
		String all_c_input_file_name__set_D = "";
		String[] all_c_input_file_name = { "russia_march.tsv" };
		//
		//
		all_c_input_file_name__set_T = c_commoon_path + "/" + "query_russia_march_1__SORTED__TOP__TSPR.tsv";
		all_c_input_file_name__set_D = c_commoon_path + "/" + "query_russia_march_2__SORTED__TOP__1_minus_TSPR.tsv";

//		all_c_input_file_name__set_T = c_commoon_path + "/" + "query_russia_march_2__SORTED__TOP__TSPR.tsv";
//		all_c_input_file_name__set_D = c_commoon_path + "/" + "query_russia_march_1__SORTED__TOP__1_minus_TSPR.tsv";

		double lambda_1 = 0.;
		double lambda_2 = 100000.;

		String c_input_file_name;
		for (String input_file_name : all_c_input_file_name) {
			c_input_file_name = c_commoon_path + "/" + input_file_name;
			//
			//
			//
			System.out.println();
			System.out.println("-----------------------------------------");
			System.out.println(c_input_file_name);
			System.out.println();
			//
//			System.gc();
			SuperGreedyPP neg_dsg_c = new SuperGreedyPP(c_input_file_name);
			neg_dsg_c.createGraph();
			neg_dsg_c.createSetsOfNodesToConsider(all_c_input_file_name__set_T, all_c_input_file_name__set_D, lambda_1,
					lambda_2);
//			NegDSGCharikarBabis neg_dsg_c = loadSerializeGraphFromFile(c_input_file_name);
//			System.gc();
			//
			//
			//
			System.out.println("input_file_name: " + input_file_name);
			//

//			int T = 10000;
			int T = 1000;
//			int T = 100;
//			int T = 10;
			double min_appx_factor = 0.999;
			neg_dsg_c.runSuperGreedyPlusPlusOnlySolutionValue(T, min_appx_factor);
		}
		return;
	}

	protected static void exp_11(String algo) throws Exception {
		// keep
//		String c_commoon_path = "/Users/ikki/DATASETS/NegDSG";
//		String[] all_c_input_file_name = { "com_amazon_ungraph__FINAL.tsv" } ;
//		String[] all_c_input_file_name = { "amz_3_colors_final.tsv" };
		//
//		c_commoon_path = "/Users/ikki/DATASETS/NegDSG/vaxnovax";
//		String[] all_c_input_file_name = { "retweet_edges.tsv" };
//		String[] all_c_input_file_name = { "karate_network.tsv" };
//		String[] all_c_input_file_name = { "polbooks__network.tsv" };
		//
//		c_commoon_path = "/Users/ikki/Downloads/ISI/";
//		String[] all_c_input_file_name = { "retweet_graph_mothersday_threshold_largest_CC.tsv" };

//		String all_c_input_file_name_NODES_TO_CONSIDER_SET_1 = c_commoon_path + "/" + "Camera_and_Photo.tsv";
//		String all_c_input_file_name_NODES_TO_CONSIDER_SET_2 = c_commoon_path + "/" + "Musical_Instruments.tsv";
		// retweet_edges.tsv
//		String all_c_input_file_name_NODES_TO_CONSIDER_SET_1 = c_commoon_path + "/"
//				+ "metis_score.csv__NO_VAX_SCORE.tsv";
//		String all_c_input_file_name_NODES_TO_CONSIDER_SET_2 = c_commoon_path + "/"
//				+ "metis_score.csv__PRO_VAX_SCORE.tsv";
//		String all_c_input_file_name_NODES_TO_CONSIDER_SET_1 = c_commoon_path + "/"
//				+ "metis_score.csv__PRO_VAX_SCORE__01.tsv";
//		String all_c_input_file_name_NODES_TO_CONSIDER_SET_2 = all_c_input_file_name_NODES_TO_CONSIDER_SET_1; 
		//
		//
//		String all_c_input_file_name_NODES_TO_CONSIDER_SET_1 = c_commoon_path + "/"
//				+ "karate_Mr_Hi.tsv";
//		String all_c_input_file_name_NODES_TO_CONSIDER_SET_2 = c_commoon_path + "/"
//				+ "karate_Officer.tsv";
		//
//		String all_c_input_file_name_NODES_TO_CONSIDER_SET_1 = c_commoon_path + "/" + "polbooks__C.tsv";
//		String all_c_input_file_name_NODES_TO_CONSIDER_SET_2 = c_commoon_path + "/" + "polbooks__not_C.tsv";
//		String all_c_input_file_name_NODES_TO_CONSIDER_SET_1 = c_commoon_path + "/" + "polbooks__L.tsv";
//		String all_c_input_file_name_NODES_TO_CONSIDER_SET_2 = c_commoon_path + "/" + "polbooks__not_L.tsv";
		//
//		String all_c_input_file_name_NODES_TO_CONSIDER_SET_1 = c_commoon_path + "/"
//				+ "retweet_graph_mothersday_threshold_largest_CC__complementingNodes__1643135238869.tsv";
//		String all_c_input_file_name_NODES_TO_CONSIDER_SET_1 = c_commoon_path + "/" + "dummy.tsv";
//		String all_c_input_file_name_NODES_TO_CONSIDER_SET_2 = c_commoon_path + "/" + "dummy.tsv";

		String c_commoon_path = "/Users/ikki/Downloads/ISI/retweet_russia_march";
		String[] all_c_input_file_name = { "retweet_graph_russia_march_threshold_largest_CC.tsv" };
		String all_c_input_file_name__set_T = "community1_russia_march__TSPR.tsv";
		String all_c_input_file_name__set_D = "community2_russia_march__1_minus_TSPR.tsv";

		double lambda_1 = 0.;
		double lambda_2 = 0.;

		String c_input_file_name;
		for (String input_file_name : all_c_input_file_name) {
			c_input_file_name = c_commoon_path + "/" + input_file_name;
			//
			//
			//
			System.out.println();
			System.out.println("-----------------------------------------");
			System.out.println(c_input_file_name);
			System.out.println();
			//
//			System.gc();
			SuperGreedyPP neg_dsg_c = new SuperGreedyPP(c_input_file_name);
			neg_dsg_c.createGraph();
			neg_dsg_c.createSetsOfNodesToConsider(all_c_input_file_name__set_T, all_c_input_file_name__set_D, lambda_1,
					lambda_2);
//			NegDSGCharikarBabis neg_dsg_c = loadSerializeGraphFromFile(c_input_file_name);
//			System.gc();
			//
			//
			//
			System.out.println("input_file_name: " + input_file_name);
			//

//			int T = 10000;
			int T = 1000;
//			int T = 100;
//			int T = 10;
			double min_appx_factor = 1; // 0.999;
//			neg_dsg_c.runSuperGreedyPlusPlusOnlySolutionValue(T, min_appx_factor);

			neg_dsg_c.lambda_1 = 0;
			neg_dsg_c.lambda_2 = 0;
			neg_dsg_c.runSuperGreedyPlusPlusOnlySolutionValue(T, min_appx_factor);
			for (int exponent = 16; exponent <= 16; exponent++) {
				int l_1 = (int) Math.pow(2, exponent);
				neg_dsg_c.lambda_1 = l_1;
				neg_dsg_c.lambda_2 = 0;
				neg_dsg_c.runSuperGreedyPlusPlusOnlySolutionValue(T, min_appx_factor);
			}
//			for (int l_1 = 0; l_1 <= 100; l_1++) {
//				for (int l_2 = 0; l_2 <= 100; l_2++) {
//					//
//					neg_dsg_c.lambda_1 = l_1;
//					neg_dsg_c.lambda_2 = l_2;
//					neg_dsg_c.runSuperGreedyPlusPlusOnlySolutionValue(T, min_appx_factor);
//					//
//				}
//			}

			//

			//
			//
			//
		}
		return;
	}

	protected static void exp_10(String algo) throws Exception {
		// keep
		String c_commoon_path = "/Users/ikki/DATASETS/NegDSG";
		String[] all_c_input_file_name = { "com_amazon_ungraph__FINAL.tsv" };
//		String[] all_c_input_file_name = { "amz_3_colors_final.tsv" };

		String c_input_file_name;
		for (String input_file_name : all_c_input_file_name) {
			c_input_file_name = c_commoon_path + "/" + input_file_name;
			//
			//
			//
			System.out.println();
			System.out.println("-----------------------------------------");
			System.out.println(c_input_file_name);
			System.out.println();
			//
//			System.gc();
			SuperGreedyPP neg_dsg_c = new SuperGreedyPP(c_input_file_name);
			neg_dsg_c.createGraph();
//			NegDSGCharikarBabis neg_dsg_c = loadSerializeGraphFromFile(c_input_file_name);
//			System.gc();
			//
			//
			//
			System.out.println("input_file_name: " + input_file_name);
			//

//			int T = 10000;
			int T = 100;
			double min_appx_factor = 0.99;
			neg_dsg_c.runSuperGreedyPlusPlusOnlySolutionValue(T, min_appx_factor);

			//

			//
			//
			//
		}
		return;
	}

	/**
	 * experiment with planted clique.
	 * 
	 * @param c_input_file_name
	 * @throws Exception
	 */
	protected static void exp_8(String c_input_file_name) throws Exception {
		//
		//
		//
		System.out.println();
		System.out.println("-----------------------------------------");
		System.out.println(c_input_file_name);
		System.out.println();
		//
		System.gc();
		// NegDSGCharikarBabis neg_dsg_c = new NegDSGCharikarBabis(c_input_file_neme);
		// neg_dsg_c.createGraph();
		SuperGreedyPP neg_dsg_c = loadSerializeGraphFromFile(c_input_file_name);
		System.gc();
		//
		//
		//
//		int T = 10000;
		int T = 1000;
		//
		//
		//
		double min_appx_factor = 0.99;
		neg_dsg_c.runSuperGreedyPlusPlusOnlySolutionValue(T, min_appx_factor);
		//
		return;
	}

	public static void serializeGraphOnFile(String input_file_complete_path) throws Exception {
		//
		SuperGreedyPP neg_dsg_c = new SuperGreedyPP(input_file_complete_path);
		neg_dsg_c.createGraph();
		//
		FileOutputStream fos = new FileOutputStream(
				input_file_complete_path.replace(".tsv", "__SERIALIZED_JAVA_OBJECT.obj"));
//		FileOutputStream fos = new FileOutputStream(
//				"/Volumes/AdriEHD4TBr/NEGDSG/DATASETS/__SERIALIZED_JAVA_OBJECT.obj");

		ObjectOutputStream oos = new ObjectOutputStream(fos);
		oos.writeObject(neg_dsg_c);
		oos.close();
		//
		return;
	}

	public static SuperGreedyPP loadSerializeGraphFromFile(String input_file_complete_path) throws Exception {
		//
		FileInputStream fis = new FileInputStream(input_file_complete_path);
		BufferedInputStream bif = new BufferedInputStream(fis, 1000000000);
		ObjectInputStream ois = new ObjectInputStream(bif);
		Object obj = ois.readObject();
		SuperGreedyPP neg_dsg_c = (SuperGreedyPP) obj;
		ois.close();
		//
		return neg_dsg_c;
	}

	// ------------------------------------------------------------------------
	// ------------------------------------------------------------------------
	// ------------------------------------------------------------------------

	public void printGraphOnScreen() {
		System.out.println();
		System.out.println("graph_file_name " + this.graph_file_name);
		System.out.println();
		for (int i = 0; i < this.graph_as_adjacency_list.length; i++) {
			System.out.println(i + " -- " + Arrays.toString(this.graph_as_adjacency_list[i]));
			System.out.println(i + " -- " + Arrays.toString(this.graph_as_adjacency_list_weights[i]));
			System.out.println();
//			for (int j = 0; j < this.graph_as_adjacency_list[i].length; j++) {
//			}
		}
		System.out.println();
		return;
	}

	protected String graph_file_name;
	//
	protected String[] map__inner_id__outer_id;
	//
	protected int[][] graph_as_adjacency_list;
	protected double[][] graph_as_adjacency_list_weights;
	protected double E_w;
	protected int V;
	protected double[] map__node__weighted_degree;
	//
	protected int max_node_id;
	protected int total_number_of_edges;
	//
	protected int output_best_V_size = -1;
	protected double output_best_solution_value = Double.NEGATIVE_INFINITY;
	protected double output_best_E_w_value = Double.NEGATIVE_INFINITY;
	protected boolean[] output_best_V = null;
	protected int index_for_best_solution_set_of_nodes = 0;
	protected int[] list__removed_nodes = null;
	protected double output_UB_as_max_inner_degree = Double.POSITIVE_INFINITY;
	protected double output_UB = Double.POSITIVE_INFINITY;
	//
	protected double[] map__node__weight_1;
	protected double[] map__node__weight_2;
	protected double sum_all_weight_1;
	protected double sum_all_weight_2;
	protected double lambda_0 = 1.;
	protected double lambda_1;
	protected double lambda_2;
	//
	int[] map__old_node_id__new_node_id = null;
	//
	protected Map<String, String> map__field__value = null;
	//
	static double DSG_density = 0.;
	//

	// ------------------------------------------------------------------------
	// ------------------------------------------------------------------------
	// ------------------------------------------------------------------------

	public void resetAllOutputAttributes() {
		this.output_best_V_size = -1;
		this.output_best_solution_value = Double.NEGATIVE_INFINITY;
		this.output_best_E_w_value = Double.NEGATIVE_INFINITY;
		this.output_best_V = null;
		return;
	}

	public int getBestVsize() {
		return this.output_best_V_size;
	}

	public double getBestDensity() {
		return this.output_best_solution_value;
	}

	public double getBestEwValue() {
		return this.output_best_E_w_value;
	}

	public boolean[] getBestV() {
		return this.output_best_V;
	}

	public SuperGreedyPP(String graph_file_name) {
		this.graph_file_name = graph_file_name;
		return;
	}

	/**
	 * 
	 * @throws Exception
	 */
	protected void createGraph() throws Exception {
		//
		Map<String, Integer> map__outer_id__inner_id = new HashMap<String, Integer>();
		String line;
		String[] record;
		BufferedReader br = new BufferedReader(new FileReader(this.graph_file_name), 100000000);
		this.total_number_of_edges = 0;
		int c_node_id = 0;
		int i;
		while ((line = br.readLine()) != null) {
			record = line.split("\t");
			//
			this.total_number_of_edges++;
			//
			for (i = 0; i < 2; i++) {
				if (!map__outer_id__inner_id.containsKey(record[i])) {
					map__outer_id__inner_id.put(record[i], c_node_id);
					c_node_id++;
				}
			}
			//
		}
		br.close();
		//
		this.max_node_id = map__outer_id__inner_id.size() - 1;
		this.map__inner_id__outer_id = new String[map__outer_id__inner_id.size()];
		Iterator<Map.Entry<String, Integer>> it = map__outer_id__inner_id.entrySet().iterator();
		Map.Entry<String, Integer> c_pair = null;
		while (it.hasNext()) {
			c_pair = it.next();
			//
			this.map__inner_id__outer_id[c_pair.getValue()] = c_pair.getKey();
			//
		}
		//
		//
		System.out.println("total #edges= " + this.total_number_of_edges);
		System.out.println("max_node_id= " + this.max_node_id);
		//
		// map__node__num_neighbours
		System.out.println("Start phase 2.");
		int c_node_id_1, c_node_id_2;
		int[] map__node__num_neighbours = new int[max_node_id + 1];
		br = new BufferedReader(new FileReader(this.graph_file_name), 100000000);
		// while ((line = csv_reader.readNext()) != null) {
		while ((line = br.readLine()) != null) {
			record = line.split("\t");
			//
			c_node_id_1 = map__outer_id__inner_id.get(record[0]);
			c_node_id_2 = map__outer_id__inner_id.get(record[1]);
			//
			map__node__num_neighbours[c_node_id_1]++;
			map__node__num_neighbours[c_node_id_2]++;
			//
			// if (++num_edges % 1000 == 0) {
			// System.out.println("line_number : " + num_edges);
			// }
		}
		br.close();
		System.out.println("end phase 2.");

		// adjacency list creation
		System.out.println("Start phase 3.");
		this.graph_as_adjacency_list = new int[this.max_node_id + 1][];
		this.graph_as_adjacency_list_weights = new double[this.max_node_id + 1][];
		this.V = 0;
		int node_id;
		for (node_id = 0; node_id < map__node__num_neighbours.length; node_id++) {
			if (map__node__num_neighbours[node_id] > 0) {
				this.V++;
				this.graph_as_adjacency_list[node_id] = new int[map__node__num_neighbours[node_id]];
				this.graph_as_adjacency_list_weights[node_id] = new double[map__node__num_neighbours[node_id]];
			}
		}
		map__node__num_neighbours = null;
		System.gc();
		//
		// this.map__node__weighted_degree = new double[this.max_node_id + 1];
		this.map__node__weighted_degree = new double[this.max_node_id + 1];
		int[] map__node__first_free_index = new int[max_node_id + 1];
		int free_index;
		double c_weight;
		br = new BufferedReader(new FileReader(this.graph_file_name), 100000000);
		this.E_w = 0;
		long num_edges = 0;
		while ((line = br.readLine()) != null) {
			record = line.split("\t");
			//
			c_node_id_1 = map__outer_id__inner_id.get(record[0]);
			c_node_id_2 = map__outer_id__inner_id.get(record[1]);
			c_weight = 1.;
			if (record.length > 2) {
				c_weight = Double.parseDouble(record[2]);
			}
			//
			this.E_w += c_weight;
			//
			// -----------------------------------------------------------------
			map__node__weighted_degree[c_node_id_1] += c_weight;
			map__node__weighted_degree[c_node_id_2] += c_weight;
			// -----------------------------------------------------------------
			//
			// this.map__node__weighted_degree[c_node_id_1] += c_weight;
			// this.map__node__weighted_degree[c_node_id_2] += c_weight;
			//
			free_index = map__node__first_free_index[c_node_id_1]++;
			this.graph_as_adjacency_list[c_node_id_1][free_index] = c_node_id_2;
			this.graph_as_adjacency_list_weights[c_node_id_1][free_index] = c_weight;
			//
			free_index = map__node__first_free_index[c_node_id_2]++;
			this.graph_as_adjacency_list[c_node_id_2][free_index] = c_node_id_1;
			this.graph_as_adjacency_list_weights[c_node_id_2][free_index] = c_weight;
			//
			if (++num_edges % 1000000 == 0) {
				System.out.println("line_number : " + num_edges);
			}
		}
		br.close();
		System.out.println("end phase 3.");

		return;
	}

	protected void createSetsOfNodesToConsider(String all_c_input_file_name_NODES_TO_CONSIDER_SET_1,
			String all_c_input_file_name_NODES_TO_CONSIDER_SET_2, double lambda_1, double lambda_2) throws Exception {
		//
		this.map__node__weight_1 = new double[this.max_node_id + 1];
		this.map__node__weight_2 = new double[this.max_node_id + 1];
		Arrays.fill(map__node__weight_1, 0.);
		Arrays.fill(map__node__weight_2, 0.);
		this.sum_all_weight_1 = 0;
		this.sum_all_weight_2 = 0;
		this.lambda_1 = lambda_1;
		this.lambda_2 = lambda_2;
		//
		HashMap<String, Integer> map__outer_id__inner_id = new HashMap<String, Integer>(2 * this.max_node_id);
		String outer_node_id = "";
		int inner_node_id = 0;
		for (inner_node_id = 0; inner_node_id < this.map__inner_id__outer_id.length; inner_node_id++) {
			outer_node_id = this.map__inner_id__outer_id[inner_node_id];
			map__outer_id__inner_id.put(outer_node_id, inner_node_id);
		}
		//
		double[][] temp_1 = new double[2][];
		temp_1[0] = this.map__node__weight_1;
		temp_1[1] = this.map__node__weight_2;
		String[] temp_2 = { all_c_input_file_name_NODES_TO_CONSIDER_SET_1,
				all_c_input_file_name_NODES_TO_CONSIDER_SET_2 };
		double[] sum_all_weights = { 0., 0. };
		//
		for (int i = 0; i < temp_1.length; i++) {
			//
			String c_file_name = temp_2[i];
			double[] map__node__weight = temp_1[i];
			//
			BufferedReader br = new BufferedReader(new FileReader(c_file_name), 1000000);
			String line = null;
			String[] record = null;
			double weight_node_id = 0.;
			while ((line = br.readLine()) != null) {
				//
				record = line.split("\t");
				outer_node_id = record[0];
				weight_node_id = Double.parseDouble(record[1]);
				//
				try {
					inner_node_id = map__outer_id__inner_id.get(outer_node_id);
				} catch (Exception e) {
					System.out.println("not alligned outer_node_id= " + outer_node_id);
					System.exit(-1);
				}
				//
				map__node__weight[inner_node_id] = weight_node_id;
				//
				sum_all_weights[i] += weight_node_id;
			}
			br.close();
		}
		//
		this.sum_all_weight_1 = sum_all_weights[0];
		this.sum_all_weight_2 = sum_all_weights[1];
		//
//		double sum = 0.;
//		for (int i = 0; i < this.map__node__weight_1.length; i++) {
//			sum += this.map__node__weight_1[i];
//		}
//		if (sum != this.sum_all_weight_1) {
//			System.out.println("sum != this.sum_all_weight_1");
//			System.out.println("sum                   "  + sum);
//			System.out.println("this.sum_all_weight_1 "  + this.sum_all_weight_1);
//			System.exit(-1);
//		}
//		sum = 0.;
//		for (int i = 0; i < this.map__node__weight_2.length; i++) {
//			sum += this.map__node__weight_2[i];
//		}
//		if (sum != this.sum_all_weight_2) {
//			System.out.println("sum != this.sum_all_weight_2");
//			System.exit(-1);
//		}
		//
		return;
	}

	protected double computeSubgraphAverageNodeAttribute(boolean[] subgraph_as_flag_array,
			double[] map__node__attribute) {
		double sum = 0.;
		double num_nodes_in_subgraph = 0.;
		for (int node_id = 0; node_id < subgraph_as_flag_array.length; node_id++) {
			if (subgraph_as_flag_array[node_id]) {
				sum += map__node__attribute[node_id];
				num_nodes_in_subgraph++;
			}
		}
		return (sum / num_nodes_in_subgraph);
	}

	protected double computeSubgraphProximityToTau(boolean[] subgraph_as_flag_array) {
		return this.computeSubgraphAverageNodeAttribute(subgraph_as_flag_array, this.map__node__weight_1);
	}

	protected double computeSubgraphDistanceFromDanger(boolean[] subgraph_as_flag_array) {
		return this.computeSubgraphAverageNodeAttribute(subgraph_as_flag_array, this.map__node__weight_2);
	}

	protected double computeSubgraphObjectiveFunctionValue(boolean[] subgraph_as_flag_array) {
		double of_value = 0.;
		of_value += this.lambda_0 * this.computeSubgraphDensity(subgraph_as_flag_array);
		of_value += this.lambda_1 * this.computeSubgraphProximityToTau(subgraph_as_flag_array);
		of_value += this.lambda_2 * this.computeSubgraphDistanceFromDanger(subgraph_as_flag_array);
		return of_value;
	}

	protected double computeSubgraphDensity(boolean[] subgraph_as_flag_array) {
		double density = 0.;
		int subgraph_num_nodes = 0;
		double subgraph_SUM_edges_weights = 0.;
		//
		double c_node_subgraph_induced_WEIGHTED_degree;
		int c_node_inner_id;
		int j;
		int c_node_NEIG_inner_id;
		int[] c_node_adj_list;
		double[] c_node_adj_list_weights;
		for (c_node_inner_id = 0; c_node_inner_id < this.graph_as_adjacency_list.length; c_node_inner_id++) {
			if (!subgraph_as_flag_array[c_node_inner_id]) {
				continue;
			}
			//
			subgraph_num_nodes++;
			//
			c_node_adj_list = this.graph_as_adjacency_list[c_node_inner_id];
			c_node_adj_list_weights = this.graph_as_adjacency_list_weights[c_node_inner_id];
			c_node_subgraph_induced_WEIGHTED_degree = 0.;
			for (j = 0; j < c_node_adj_list.length; j++) {
				c_node_NEIG_inner_id = c_node_adj_list[j];
				if (!subgraph_as_flag_array[c_node_NEIG_inner_id]) {
					continue;
				}
				//
				c_node_subgraph_induced_WEIGHTED_degree += c_node_adj_list_weights[j];
			}
			//
			subgraph_SUM_edges_weights += c_node_subgraph_induced_WEIGHTED_degree;
			//
		}
		//
		subgraph_SUM_edges_weights /= 2;
		density = subgraph_SUM_edges_weights / subgraph_num_nodes;
		//
		return density;
	}

	protected int computeSubgraphNumNodes(boolean[] subgraph_as_flag_array) {
		int num_nodes = 0;
		for (int i = 0; i < subgraph_as_flag_array.length; i++) {
			if (subgraph_as_flag_array[i]) {
				num_nodes++;
			}
		}
		return num_nodes;
	}

	protected double computeSumAllEdgesWeights() {
		double graph_SUM_edges_weights = 0.;
		//
		double c_node_subgraph_induced_WEIGHTED_degree;
		int c_node_inner_id;
		int j;
		int[] c_node_adj_list;
		double[] c_node_adj_list_weights;
		for (c_node_inner_id = 0; c_node_inner_id < this.graph_as_adjacency_list.length; c_node_inner_id++) {
			//
			c_node_adj_list = this.graph_as_adjacency_list[c_node_inner_id];
			c_node_adj_list_weights = this.graph_as_adjacency_list_weights[c_node_inner_id];
			c_node_subgraph_induced_WEIGHTED_degree = 0.;
			for (j = 0; j < c_node_adj_list.length; j++) {
				c_node_subgraph_induced_WEIGHTED_degree += c_node_adj_list_weights[j];
			}
			//
			graph_SUM_edges_weights += c_node_subgraph_induced_WEIGHTED_degree;
			//
		}
		//
		graph_SUM_edges_weights /= 2.;
		//
		return graph_SUM_edges_weights;
	}

	protected int computeSubgraphNumberOfNodes(boolean[] subgraph_as_flag_array) {
		int num_nodes = 0;
		//
		for (int i = 0; i < subgraph_as_flag_array.length; i++) {
			if (subgraph_as_flag_array[i]) {
				num_nodes++;
			}
		}
		//
		return num_nodes;
	}

	protected String convertSubgraphInSring(boolean[] subgraph_as_flag_array) {
		String out_string = "";
		//
		StringBuilder builder = new StringBuilder();
		boolean fist_node = true;
		for (int i = 0; i < subgraph_as_flag_array.length; i++) {
			if (subgraph_as_flag_array[i]) {
				if (fist_node) {
					fist_node = false;
					builder.append("{");
				} else {
					builder.append(", ");
				}
				builder.append(this.map__inner_id__outer_id[i]);
			}
		}
		builder.append("}");
		out_string = builder.toString();
		//
		return out_string;
	}

	protected String convertSubgraphInSring2(boolean[] subgraph_as_flag_array) {
		String out_string = "";
		//
		ArrayList<Integer> lllll = new ArrayList<Integer>();
		int e;
		for (int i = 0; i < subgraph_as_flag_array.length; i++) {
			if (subgraph_as_flag_array[i]) {
				e = Integer.parseInt(this.map__inner_id__outer_id[i]);
				lllll.add(e);
			}
		}
		Collections.sort(lllll);
		out_string = lllll.toString();
		//
		return out_string;
	}

	protected String convertSubgraphInSring3(boolean[] subgraph_as_flag_array) {
		String out_string = "";
		//
		ArrayList<String> lllll = new ArrayList<String>();
		String e;
		for (int i = 0; i < subgraph_as_flag_array.length; i++) {
			if (subgraph_as_flag_array[i]) {
				e = this.map__inner_id__outer_id[i];
				lllll.add(e);
			}
		}
		Collections.sort(lllll);
		out_string = lllll.toString();
		// out_string = out_string.replace("[", "{").replace("]", "}");
		//
		return out_string;
	}

	protected String printComplementingNodes(boolean[] subgraph_as_flag_array) {
		String out_string = "";
		//
		StringBuilder builder = new StringBuilder();
		for (int i = 0; i < subgraph_as_flag_array.length; i++) {
			if (!subgraph_as_flag_array[i]) {
				builder.append(this.map__inner_id__outer_id[i]);
				builder.append("\t1\n");
			}
		}
		out_string = builder.toString();
		//
		return out_string;
	}

	protected double[] computeSubgraphBalance(boolean[] subgraph_as_flag_array) {
		double[] map__class__num_nodes = { 0., 0., 0. };
		//
		for (int i = 0; i < subgraph_as_flag_array.length; i++) {
			if (subgraph_as_flag_array[i]) {
				map__class__num_nodes[0] += this.map__node__weight_1[i];
				map__class__num_nodes[1] += this.map__node__weight_2[i];
				if (this.map__node__weight_1[i] + this.map__node__weight_2[i] == 0.) {
					map__class__num_nodes[2] += 1;
				}
			}
		}
		//
		return map__class__num_nodes;
	}

	protected SuperGreedyPP createSubGraph(boolean[] map__node_id__deleted_node) throws Exception {
		//
//		System.out.println();
//		System.out.println(Arrays.toString(map__node_id__deleted_node));
//		System.out.println();
		//
		SuperGreedyPP subg = new SuperGreedyPP(this.graph_file_name + "__" + "subgraph");
		//
		int j, j_sub, c_node_id, prev_node_id;
		int c_node_id_i, c_node_id_j, prev_node_id_i, prev_node_id_j;
		double c_node_id_j__weight;
		//
		this.map__old_node_id__new_node_id = new int[this.max_node_id + 1];
		c_node_id = -1;
		for (prev_node_id = 0; prev_node_id < map__node_id__deleted_node.length; prev_node_id++) {
			map__old_node_id__new_node_id[prev_node_id] = -1;
			if (!map__node_id__deleted_node[prev_node_id]) {
				c_node_id++;
				map__old_node_id__new_node_id[prev_node_id] = c_node_id;
			}
		}
		subg.max_node_id = c_node_id;
		subg.V = subg.max_node_id + 1;
		//
		// System.out.println("subg.max_node_id " + subg.max_node_id);
		//
		subg.map__inner_id__outer_id = new String[subg.max_node_id + 1];

		c_node_id = -1;
		for (prev_node_id = 0; prev_node_id < this.map__inner_id__outer_id.length; prev_node_id++) {
			//
			c_node_id = map__old_node_id__new_node_id[prev_node_id];
			if (c_node_id < 0) {
				continue;
			}
			//
			// System.out.println("c_node_id: " + c_node_id);
			//
			map__old_node_id__new_node_id[prev_node_id] = c_node_id;
			//
			subg.map__inner_id__outer_id[c_node_id] = this.map__inner_id__outer_id[prev_node_id];
			//
		}
		//
		//
		subg.graph_as_adjacency_list = new int[subg.max_node_id + 1][subg.max_node_id + 1];
		subg.graph_as_adjacency_list_weights = new double[subg.max_node_id + 1][subg.max_node_id + 1];
		subg.map__node__weighted_degree = new double[subg.max_node_id + 1];
		subg.E_w = 0.;
		subg.total_number_of_edges = 0;
		for (prev_node_id_i = 0; prev_node_id_i < this.graph_as_adjacency_list.length; prev_node_id_i++) {
			c_node_id_i = map__old_node_id__new_node_id[prev_node_id_i];
			if (c_node_id_i < 0) {
				continue;
			}
			//
			j_sub = -1;
			for (j = 0; j < this.graph_as_adjacency_list[prev_node_id_i].length; j++) {
				prev_node_id_j = this.graph_as_adjacency_list[prev_node_id_i][j];
				c_node_id_j = map__old_node_id__new_node_id[prev_node_id_j];
				if (c_node_id_j < 0) {
					continue;
				}
				j_sub++;
				//
				c_node_id_j__weight = this.graph_as_adjacency_list_weights[prev_node_id_i][j];
				subg.E_w += c_node_id_j__weight;
				subg.total_number_of_edges++;
				//
				subg.graph_as_adjacency_list[c_node_id_i][j_sub] = c_node_id_j;
				subg.graph_as_adjacency_list_weights[c_node_id_i][j_sub] = c_node_id_j__weight;
				//
				subg.map__node__weighted_degree[c_node_id_i] += c_node_id_j__weight;
				//
			}
		}
		subg.E_w /= 2.;
		subg.total_number_of_edges /= 2;
		subg.output_best_V_size = -1;
		subg.output_best_solution_value = Double.NEGATIVE_INFINITY;
		subg.output_best_E_w_value = Double.NEGATIVE_INFINITY;
		subg.output_best_V = null;
		//
		subg.map__node__weight_1 = new double[subg.max_node_id + 1];
		subg.map__node__weight_2 = new double[subg.max_node_id + 1];
		subg.lambda_0 = 1.;
		subg.lambda_1 = 0.;
		subg.lambda_2 = 0.;
		//
		return subg;
	}

	public boolean[] BACKUP__singleIterationOfSuperGreedyPlusPlusOnlySolutionValue(double[] map__node_id__load,
			boolean complete_output) {
		// keep
		int c_node_id;
		//
		boolean[] map__node__deleted = new boolean[this.max_node_id + 1];
		// Arrays.fill(map__node__deleted, false);
		//
		double c_load = 0.;
		double min_load = Double.POSITIVE_INFINITY;
		double max_load = Double.NEGATIVE_INFINITY;
		double c_inner_degree = 0.;
		double max_inner_degree = Double.NEGATIVE_INFINITY;
		//
		int index_for_best_solution_set_of_nodes = 0;
		int[] list__removed_nodes = null;
		int num_nodes_removed_so_far = 0;
		if (complete_output) {
			list__removed_nodes = new int[this.max_node_id + 1];
			Arrays.fill(list__removed_nodes, -1);
		}
		//
		// "PQ"
		double c_node_priority;
		FastAndSimpleIndexMinPQ pq = new FastAndSimpleIndexMinPQ(this.map__node__weighted_degree.length);
		for (c_node_id = 0; c_node_id < this.map__node__weighted_degree.length; c_node_id++) {
			//
			c_node_priority = map__node_id__load[c_node_id];
			c_node_priority += this.map__node__weighted_degree[c_node_id]; // ric
			//
			pq.insert(c_node_id, c_node_priority);
		}
		//
		double c_E_w_value = this.E_w;
		int c_V_value = this.V;
		double c_density = Double.NEGATIVE_INFINITY;
		//
		int c_highest_priority_node_id = -1;
		int i, c_not_removed_neighbour_id;
		int[] list_of_all_original_neighbours_of_c_highest_priority_node_id;
		double[] list_of_weights_of_all_original_neighbours_of_c_highest_priority_node_id;
		double PRIORITY_of_c_not_removed_neighbour_of_c_min_weighted_degree_node;
		double new_PRIORITY_of_c_not_removed_neighbour_of_c_min_weighted_degree_node;
		double weight_on_current_edge;
		//
		double best_solution__density = -1;
		int best_solution__size = -1;
		while (c_V_value >= 2) {
			//
			// compute density of current induced subgraph
			c_density = ((double) (c_E_w_value)) / c_V_value;
			//
			if (c_density >= best_solution__density) {
				best_solution__density = c_density;
				best_solution__size = c_V_value;
				if (complete_output) {
					index_for_best_solution_set_of_nodes = num_nodes_removed_so_far;
				}
			}
			//
			// save the current configuration
			if (complete_output) {
//				list__removed_nodes_density[num_nodes_removed_so_far] = c_density;
				list__removed_nodes[num_nodes_removed_so_far] = c_highest_priority_node_id;
				num_nodes_removed_so_far++;
			}
			//
			// Compute the new load for the node with minimum weighted degree in the
			// priority queue.
			c_load = pq.minKey();
			//
			min_load = (c_load < min_load ? c_load : min_load);
			max_load = (c_load > max_load ? c_load : max_load);
			//
			// extract the node with minimum weighted degree in the priority queue.
			c_highest_priority_node_id = pq.delMin();
			//
			// Update the load for the node with minimum weighted degree in the
			// priority queue.
			c_inner_degree = c_load - map__node_id__load[c_highest_priority_node_id];
			max_inner_degree = (c_inner_degree > max_inner_degree ? c_inner_degree : max_inner_degree);
			//
			// remove the node with minimum weighted degree from the graph
			map__node__deleted[c_highest_priority_node_id] = true;
			//
			// update priorities on priority queue
			// iterate over all not removed neighbours of
			// c_highest_priority_node.
			list_of_all_original_neighbours_of_c_highest_priority_node_id = this.graph_as_adjacency_list[c_highest_priority_node_id];
			list_of_weights_of_all_original_neighbours_of_c_highest_priority_node_id = this.graph_as_adjacency_list_weights[c_highest_priority_node_id];
			for (i = 0; i < list_of_all_original_neighbours_of_c_highest_priority_node_id.length; i++) {
				c_not_removed_neighbour_id = list_of_all_original_neighbours_of_c_highest_priority_node_id[i];
				if (map__node__deleted[c_not_removed_neighbour_id]) {
					continue;
				}
				//
				// Update priority of the not_removed_neighbour_of_c_min_weighted_degree_node
				PRIORITY_of_c_not_removed_neighbour_of_c_min_weighted_degree_node = pq
						.keyOf(c_not_removed_neighbour_id);
				weight_on_current_edge = list_of_weights_of_all_original_neighbours_of_c_highest_priority_node_id[i];
				if (weight_on_current_edge >= 0) {
					new_PRIORITY_of_c_not_removed_neighbour_of_c_min_weighted_degree_node = PRIORITY_of_c_not_removed_neighbour_of_c_min_weighted_degree_node
							- weight_on_current_edge;
					//
					map__node_id__load[c_highest_priority_node_id] += weight_on_current_edge;
					//
				} else {
					new_PRIORITY_of_c_not_removed_neighbour_of_c_min_weighted_degree_node = PRIORITY_of_c_not_removed_neighbour_of_c_min_weighted_degree_node
							- weight_on_current_edge;
				}
				//
				pq.updateKey(c_not_removed_neighbour_id,
						new_PRIORITY_of_c_not_removed_neighbour_of_c_min_weighted_degree_node);
				//
				// update c_E_w_value
				c_E_w_value -= weight_on_current_edge;
			}
			//
			// update c_V_value
			c_V_value--;
		}
		//
		//
		//
		//
		if (complete_output) {
			if ((this.output_best_V == null) || (this.output_best_V.length != (this.max_node_id + 1))) {
				this.output_best_V = new boolean[this.max_node_id + 1];
			}
			Arrays.fill(this.output_best_V, true);
			for (i = 1; i <= index_for_best_solution_set_of_nodes; i++) {
				c_node_id = list__removed_nodes[i];
				this.output_best_V[c_node_id] = false;
			}
		} else {
			this.output_best_V = null;
		}
		//
		this.output_best_V_size = best_solution__size;
		this.output_best_solution_value = best_solution__density;
		this.output_best_E_w_value = best_solution__density * best_solution__size;
		//
		this.output_UB_as_max_inner_degree = max_inner_degree;
		this.output_UB = max_load / min_load;
		//
		return this.output_best_V;
	}

	public boolean[] SLOW_COMPUTATION_OF_COMPLETE_OUTPUT_singleIterationOfSuperGreedyPlusPlusOnlySolutionValue(
			double[] map__node_id__load, boolean complete_output) {
		// keep
		System.out.println("rwetrrwgtrwgsradsar");
		int c_node_id;
		//
		boolean[] map__node__deleted = new boolean[this.max_node_id + 1];
		//
//		int index_for_best_solution_set_of_nodes = 0;
//		int[] list__removed_nodes = null;
		if (complete_output) {
			this.output_best_V = new boolean[this.max_node_id + 1];
			Arrays.fill(this.output_best_V, false);
		}
		//
		// "PQ"
		double c_node_priority;
		double c_node_lambda_weighted_sum_on_it = 0.;
		FastAndSimpleIndexMinPQ pq = new FastAndSimpleIndexMinPQ(this.map__node__weighted_degree.length);
		for (c_node_id = 0; c_node_id < this.map__node__weighted_degree.length; c_node_id++) {
			//
			c_node_priority = map__node_id__load[c_node_id];
			c_node_priority += this.map__node__weighted_degree[c_node_id];
			c_node_lambda_weighted_sum_on_it = this.lambda_1 * this.map__node__weight_1[c_node_id]
					+ this.lambda_2 * this.map__node__weight_2[c_node_id];
			c_node_priority += c_node_lambda_weighted_sum_on_it;
			//
			pq.insert(c_node_id, c_node_priority);
		}
		//
		double c_E_w_value = this.E_w;
		double c_sum_weights = this.lambda_1 * this.sum_all_weight_1 + this.lambda_2 * this.sum_all_weight_2;
		int c_V_value = this.V;
		double c_solution_value = Double.NEGATIVE_INFINITY;
		int c_highest_priority_node_id = -1;
		int i, c_not_removed_neighbour_id;
		int[] list_of_all_original_neighbours_of_c_highest_priority_node_id;
		double[] list_of_weights_of_all_original_neighbours_of_c_highest_priority_node_id;
		double PRIORITY_of_c_not_removed_neighbour_of_c_min_weighted_degree_node;
		double new_PRIORITY_of_c_not_removed_neighbour_of_c_min_weighted_degree_node;
		double weight_on_current_edge;
		//
		double best_solution__value = -1;
		int best_solution__size = -1;
		while (c_V_value >= 1) {
			//
			// compute the value of current induced subgraph
			c_solution_value = (c_E_w_value + c_sum_weights) / c_V_value;
			//
			if (c_solution_value >= best_solution__value) {
				best_solution__value = c_solution_value;
				best_solution__size = c_V_value;
				if (complete_output) {
					int size__qqq = 0;
					for (int qqq = 0; qqq < map__node__deleted.length; qqq++) {
						this.output_best_V[qqq] = !map__node__deleted[qqq];
						if (this.output_best_V[qqq]) {
							size__qqq++;
						}
					}
					double density_qqq = this.computeSubgraphDensity(this.output_best_V);
					System.out.println("c_solution_size__qqq   :=" + size__qqq);
					System.out.println("best_solution__size    :=" + best_solution__size);
					System.out.println("c_solution_value__qqq  :=" + c_solution_value);
					System.out.println("c_solution_density__qqq:=" + density_qqq);
				}
			}
			//
			c_highest_priority_node_id = pq.delMin();
			//
			// Update the load for the node with minimum weighted degree in the
			// priority queue.
			c_node_lambda_weighted_sum_on_it = this.lambda_1 * this.map__node__weight_1[c_highest_priority_node_id]
					+ this.lambda_2 * this.map__node__weight_2[c_highest_priority_node_id];
			map__node_id__load[c_highest_priority_node_id] += c_node_lambda_weighted_sum_on_it;
			c_sum_weights -= c_node_lambda_weighted_sum_on_it;
			//
			// remove the node with minimum weighted degree from the graph
			map__node__deleted[c_highest_priority_node_id] = true;
			//
			// update priorities on priority queue
			// iterate over all not removed neighbours of
			// c_highest_priority_node.
			list_of_all_original_neighbours_of_c_highest_priority_node_id = this.graph_as_adjacency_list[c_highest_priority_node_id];
			list_of_weights_of_all_original_neighbours_of_c_highest_priority_node_id = this.graph_as_adjacency_list_weights[c_highest_priority_node_id];
			for (i = 0; i < list_of_all_original_neighbours_of_c_highest_priority_node_id.length; i++) {
				c_not_removed_neighbour_id = list_of_all_original_neighbours_of_c_highest_priority_node_id[i];
				if (map__node__deleted[c_not_removed_neighbour_id]) {
					continue;
				}
				//
				// Update priority of the not_removed_neighbour_of_c_min_weighted_degree_node
				PRIORITY_of_c_not_removed_neighbour_of_c_min_weighted_degree_node = pq
						.keyOf(c_not_removed_neighbour_id);
				weight_on_current_edge = list_of_weights_of_all_original_neighbours_of_c_highest_priority_node_id[i];
				new_PRIORITY_of_c_not_removed_neighbour_of_c_min_weighted_degree_node = PRIORITY_of_c_not_removed_neighbour_of_c_min_weighted_degree_node
						- weight_on_current_edge;
				//
				map__node_id__load[c_highest_priority_node_id] += weight_on_current_edge;
				//
				pq.updateKey(c_not_removed_neighbour_id,
						new_PRIORITY_of_c_not_removed_neighbour_of_c_min_weighted_degree_node);
				//
				// update c_E_w_value
				c_E_w_value -= weight_on_current_edge;
			}
			//
			// update c_V_value
			c_V_value--;
		}
		//
		if (!complete_output) {
			this.output_best_V = null;
		}
		//
		this.output_best_V_size = best_solution__size;
		this.output_best_solution_value = best_solution__value;
		//
		System.out.println("best_solution__size    |=" + best_solution__size);
		System.out.println("this.output_best_V_size|=" + this.output_best_V_size);
		//
		return this.output_best_V;
	}

	/**
	 * 
	 * @param map__node_id__load
	 * @param complete_output
	 * @return
	 */
	public boolean[] REAL__singleIterationOfSuperGreedyPlusPlusOnlySolutionValue(double[] map__node_id__load,
			boolean complete_output) {
		// keep
		int c_node_id;
		//
		boolean[] map__node__deleted = new boolean[this.max_node_id + 1];
		//
		int index_for_best_solution_set_of_nodes = 0;
		int[] list__removed_nodes = null;
		int num_nodes_removed_so_far = 0;
		if (complete_output) {
			list__removed_nodes = new int[this.max_node_id + 1];
			Arrays.fill(list__removed_nodes, -1);
		}
		//
		// "PQ"
		double c_node_priority;
		double c_node_lambda_weighted_sum_on_it = 0.;
		FastAndSimpleIndexMinPQ pq = new FastAndSimpleIndexMinPQ(this.map__node__weighted_degree.length);
		for (c_node_id = 0; c_node_id < this.map__node__weighted_degree.length; c_node_id++) {
			//
			c_node_priority = map__node_id__load[c_node_id];
			c_node_priority += this.map__node__weighted_degree[c_node_id];
			c_node_lambda_weighted_sum_on_it = this.lambda_1 * this.map__node__weight_1[c_node_id]
					+ this.lambda_2 * this.map__node__weight_2[c_node_id];
			c_node_priority += c_node_lambda_weighted_sum_on_it;
			//
			pq.insert(c_node_id, c_node_priority);
		}
		//
		double c_E_w_value = this.E_w;
		double c_sum_weights = this.lambda_1 * this.sum_all_weight_1 + this.lambda_2 * this.sum_all_weight_2;
		int c_V_value = this.V;
		double c_solution_value = Double.NEGATIVE_INFINITY;
//		System.out.println();
//		System.out.println("this.lambda_1 * this.sum_all_weight_1 " + (this.lambda_1 * this.sum_all_weight_1 ));
//		System.out.println("this.lambda_2 * this.sum_all_weight_2 "+ (this.lambda_2 * this.sum_all_weight_2));
//		System.out.println("c_sum_weights "+c_sum_weights);
		//
		int c_highest_priority_node_id = -1;
		int i, c_not_removed_neighbour_id;
		int[] list_of_all_original_neighbours_of_c_highest_priority_node_id;
		double[] list_of_weights_of_all_original_neighbours_of_c_highest_priority_node_id;
		double PRIORITY_of_c_not_removed_neighbour_of_c_min_weighted_degree_node;
		double new_PRIORITY_of_c_not_removed_neighbour_of_c_min_weighted_degree_node;
		double weight_on_current_edge;
		//
		double best_solution__value = -1;
		int best_solution__size = -1;
		while (c_V_value >= 1) {
			//
			// compute the value of current induced subgraph
			c_solution_value = (c_E_w_value + c_sum_weights) / c_V_value;
			//
			if (c_solution_value >= best_solution__value) {
				best_solution__value = c_solution_value;
				best_solution__size = c_V_value;
				if (complete_output) {
					index_for_best_solution_set_of_nodes = num_nodes_removed_so_far;
				}
			}
			//
			// save the current configuration
			if (complete_output) {
				list__removed_nodes[num_nodes_removed_so_far] = c_highest_priority_node_id;
				num_nodes_removed_so_far++;
			}
			//
			c_highest_priority_node_id = pq.delMin();
			//
			// Update the load for the node with minimum weighted degree in the
			// priority queue.
			c_node_lambda_weighted_sum_on_it = this.lambda_1 * this.map__node__weight_1[c_highest_priority_node_id]
					+ this.lambda_2 * this.map__node__weight_2[c_highest_priority_node_id];
			map__node_id__load[c_highest_priority_node_id] += c_node_lambda_weighted_sum_on_it;
			c_sum_weights -= c_node_lambda_weighted_sum_on_it;
			//
			// remove the node with minimum weighted degree from the graph
			map__node__deleted[c_highest_priority_node_id] = true;
			//
			// update priorities on priority queue
			// iterate over all not removed neighbours of
			// c_highest_priority_node.
			list_of_all_original_neighbours_of_c_highest_priority_node_id = this.graph_as_adjacency_list[c_highest_priority_node_id];
			list_of_weights_of_all_original_neighbours_of_c_highest_priority_node_id = this.graph_as_adjacency_list_weights[c_highest_priority_node_id];
			for (i = 0; i < list_of_all_original_neighbours_of_c_highest_priority_node_id.length; i++) {
				c_not_removed_neighbour_id = list_of_all_original_neighbours_of_c_highest_priority_node_id[i];
				if (map__node__deleted[c_not_removed_neighbour_id]) {
					continue;
				}
				//
				// Update priority of the not_removed_neighbour_of_c_min_weighted_degree_node
				PRIORITY_of_c_not_removed_neighbour_of_c_min_weighted_degree_node = pq
						.keyOf(c_not_removed_neighbour_id);
				weight_on_current_edge = list_of_weights_of_all_original_neighbours_of_c_highest_priority_node_id[i];
				new_PRIORITY_of_c_not_removed_neighbour_of_c_min_weighted_degree_node = PRIORITY_of_c_not_removed_neighbour_of_c_min_weighted_degree_node
						- weight_on_current_edge;
				//
				map__node_id__load[c_highest_priority_node_id] += weight_on_current_edge;
				//
				pq.updateKey(c_not_removed_neighbour_id,
						new_PRIORITY_of_c_not_removed_neighbour_of_c_min_weighted_degree_node);
				//
				// update c_E_w_value
				c_E_w_value -= weight_on_current_edge;
			}
			//
			// update c_V_value
			c_V_value--;
		}
		//
		if (complete_output) {
			if ((this.output_best_V == null) || (this.output_best_V.length != (this.max_node_id + 1))) {
				this.output_best_V = new boolean[this.max_node_id + 1];
			}
			Arrays.fill(this.output_best_V, true);
			for (i = 1; i <= index_for_best_solution_set_of_nodes; i++) {
				c_node_id = list__removed_nodes[i];
				this.output_best_V[c_node_id] = false;
			}
		} else {
			this.output_best_V = null;
		}
		//
		this.output_best_V_size = best_solution__size;
		this.output_best_solution_value = best_solution__value;
		//
		return this.output_best_V;
	}

	/**
	 * 
	 * @param map__node_id__load
	 * @param complete_output
	 * @return
	 */
	public void singleIterationOfSuperGreedyPlusPlusOnlySolutionValue(double[] map__node_id__load) {
		// keep
		int c_node_id;
		//
		boolean[] map__node__deleted = new boolean[this.max_node_id + 1];
		//
		int num_nodes_removed_so_far = 0;
		this.index_for_best_solution_set_of_nodes = 0;
		if ((this.list__removed_nodes == null) || (this.list__removed_nodes.length != (this.max_node_id + 1))) {
			this.list__removed_nodes = new int[this.max_node_id + 1];
			Arrays.fill(this.list__removed_nodes, -1);
		}
		//
		// "PQ"
		double c_node_priority;
		double c_node_lambda_weighted_sum_on_it = 0.;
		FastAndSimpleIndexMinPQ pq = new FastAndSimpleIndexMinPQ(this.map__node__weighted_degree.length);
		for (c_node_id = 0; c_node_id < this.map__node__weighted_degree.length; c_node_id++) {
			//
			c_node_priority = map__node_id__load[c_node_id];
			c_node_priority += this.lambda_0 * this.map__node__weighted_degree[c_node_id];
			c_node_lambda_weighted_sum_on_it = this.lambda_1 * this.map__node__weight_1[c_node_id]
					+ this.lambda_2 * this.map__node__weight_2[c_node_id];
			c_node_priority += c_node_lambda_weighted_sum_on_it;
			//
			pq.insert(c_node_id, c_node_priority);
		}
		//
		double c_E_w_value = this.lambda_0 * this.E_w;
		double c_sum_weights = this.lambda_1 * this.sum_all_weight_1 + this.lambda_2 * this.sum_all_weight_2;
		int c_V_value = this.V;
		double c_solution_value = Double.NEGATIVE_INFINITY;
		//
		int c_highest_priority_node_id = -1;
		int i, c_not_removed_neighbour_id;
		int[] list_of_all_original_neighbours_of_c_highest_priority_node_id;
		double[] list_of_weights_of_all_original_neighbours_of_c_highest_priority_node_id;
		double PRIORITY_of_c_not_removed_neighbour_of_c_min_weighted_degree_node;
		double new_PRIORITY_of_c_not_removed_neighbour_of_c_min_weighted_degree_node;
		double weight_on_current_edge;
		//
		double best_solution__value = -1;
		int best_solution__size = -1;
		while (c_V_value >= 1) {
			//
			// compute the value of current induced subgraph
			c_solution_value = (c_E_w_value + c_sum_weights) / c_V_value;
			//
			if (c_solution_value >= best_solution__value) {
				best_solution__value = c_solution_value;
				best_solution__size = c_V_value;
				this.index_for_best_solution_set_of_nodes = num_nodes_removed_so_far;
			}
			//
			// save the current configuration
			this.list__removed_nodes[num_nodes_removed_so_far] = c_highest_priority_node_id;
			num_nodes_removed_so_far++;

			//
			c_highest_priority_node_id = pq.delMin();
			//
			// Update the load for the node with minimum weighted degree in the
			// priority queue.
			c_node_lambda_weighted_sum_on_it = this.lambda_1 * this.map__node__weight_1[c_highest_priority_node_id]
					+ this.lambda_2 * this.map__node__weight_2[c_highest_priority_node_id];
			map__node_id__load[c_highest_priority_node_id] += c_node_lambda_weighted_sum_on_it;
			c_sum_weights -= c_node_lambda_weighted_sum_on_it;
			//
			// remove the node with minimum weighted degree from the graph
			map__node__deleted[c_highest_priority_node_id] = true;
			//
			// update priorities on priority queue
			// iterate over all not removed neighbours of
			// c_highest_priority_node.
			list_of_all_original_neighbours_of_c_highest_priority_node_id = this.graph_as_adjacency_list[c_highest_priority_node_id];
			list_of_weights_of_all_original_neighbours_of_c_highest_priority_node_id = this.graph_as_adjacency_list_weights[c_highest_priority_node_id];
			for (i = 0; i < list_of_all_original_neighbours_of_c_highest_priority_node_id.length; i++) {
				c_not_removed_neighbour_id = list_of_all_original_neighbours_of_c_highest_priority_node_id[i];
				if (map__node__deleted[c_not_removed_neighbour_id]) {
					continue;
				}
				//
				// Update priority of the not_removed_neighbour_of_c_min_weighted_degree_node
				PRIORITY_of_c_not_removed_neighbour_of_c_min_weighted_degree_node = pq
						.keyOf(c_not_removed_neighbour_id);
				weight_on_current_edge = list_of_weights_of_all_original_neighbours_of_c_highest_priority_node_id[i];
				new_PRIORITY_of_c_not_removed_neighbour_of_c_min_weighted_degree_node = PRIORITY_of_c_not_removed_neighbour_of_c_min_weighted_degree_node
						- weight_on_current_edge;
				//
				map__node_id__load[c_highest_priority_node_id] += this.lambda_0 * weight_on_current_edge;
				//
				pq.updateKey(c_not_removed_neighbour_id,
						new_PRIORITY_of_c_not_removed_neighbour_of_c_min_weighted_degree_node);
				//
				// update c_E_w_value
				c_E_w_value -= this.lambda_0 * weight_on_current_edge;
			}
			//
			// update c_V_value
			c_V_value--;
		}
		//
		//
		this.output_best_V_size = best_solution__size;
		this.output_best_solution_value = best_solution__value;
		//
		return;
	}

	protected void computeSolutionRepresentation(int index_for_best_solution_set_of_nodes, int[] list__removed_nodes) {
		if ((this.output_best_V == null) || (this.output_best_V.length != (this.max_node_id + 1))) {
			this.output_best_V = new boolean[this.max_node_id + 1];
		}
		Arrays.fill(this.output_best_V, true);
		int num_deleted_nodes = 0;
		int c_node_id;
		for (int i = 1; i <= index_for_best_solution_set_of_nodes; i++) {
			c_node_id = list__removed_nodes[i];
			this.output_best_V[c_node_id] = false;
			num_deleted_nodes++;
		}
		this.output_best_V_size = this.output_best_V.length - num_deleted_nodes;
		return;
	}

	public void singleIterationOfSuperGreedyPlusPlusOnlySolutionValue(double[] map__node_id__load,
			boolean[] map__node__tau, int num_nodes_in_tau) {
		// keep
		int c_node_id;
		//
		boolean[] map__node__deleted = new boolean[this.max_node_id + 1];
		//
		int num_nodes_removed_so_far = 0;
		this.index_for_best_solution_set_of_nodes = 0;
		if ((this.list__removed_nodes == null) || (this.list__removed_nodes.length != (this.max_node_id + 1))) {
			this.list__removed_nodes = new int[this.max_node_id + 1];
			Arrays.fill(this.list__removed_nodes, -1);
		}
		//
		// "PQ"
		double c_node_priority;
		double c_node_lambda_weighted_sum_on_it = 0.;
		FastAndSimpleIndexMinPQ pq = new FastAndSimpleIndexMinPQ(this.map__node__weighted_degree.length);
		for (c_node_id = 0; c_node_id < this.map__node__weighted_degree.length; c_node_id++) {
			//
			if (map__node__tau[c_node_id]) {
				continue;
			}
			//
			c_node_priority = map__node_id__load[c_node_id];
			c_node_priority += this.map__node__weighted_degree[c_node_id];
			c_node_lambda_weighted_sum_on_it = this.lambda_1 * this.map__node__weight_1[c_node_id]
					+ this.lambda_2 * this.map__node__weight_2[c_node_id];
			c_node_priority += c_node_lambda_weighted_sum_on_it;
			//
			pq.insert(c_node_id, c_node_priority);
		}
		//
		double c_E_w_value = this.E_w;
		double c_sum_weights = this.lambda_1 * this.sum_all_weight_1 + this.lambda_2 * this.sum_all_weight_2;
		int c_V_value = this.V;
		double c_solution_value = Double.NEGATIVE_INFINITY;
		//
		int c_highest_priority_node_id = -1;
		int i, c_not_removed_neighbour_id;
		int[] list_of_all_original_neighbours_of_c_highest_priority_node_id;
		double[] list_of_weights_of_all_original_neighbours_of_c_highest_priority_node_id;
		double PRIORITY_of_c_not_removed_neighbour_of_c_min_weighted_degree_node;
		double new_PRIORITY_of_c_not_removed_neighbour_of_c_min_weighted_degree_node;
		double weight_on_current_edge;
		//
		double best_solution__value = -1;
		int best_solution__size = -1;
		// while (c_V_value >= 1) {
		while ((c_V_value - num_nodes_in_tau) >= 1) {
			//
			// compute the value of current induced subgraph
			c_solution_value = (c_E_w_value + c_sum_weights) / c_V_value;
			//
			if (c_solution_value >= best_solution__value) {
				best_solution__value = c_solution_value;
				best_solution__size = c_V_value;
				this.index_for_best_solution_set_of_nodes = num_nodes_removed_so_far;
			}
			//
			// save the current configuration
			this.list__removed_nodes[num_nodes_removed_so_far] = c_highest_priority_node_id;
			num_nodes_removed_so_far++;

			//
			c_highest_priority_node_id = pq.delMin();
			//
			// Update the load for the node with minimum weighted degree in the
			// priority queue.
			c_node_lambda_weighted_sum_on_it = this.lambda_1 * this.map__node__weight_1[c_highest_priority_node_id]
					+ this.lambda_2 * this.map__node__weight_2[c_highest_priority_node_id];
			map__node_id__load[c_highest_priority_node_id] += c_node_lambda_weighted_sum_on_it;
			c_sum_weights -= c_node_lambda_weighted_sum_on_it;
			//
			// remove the node with minimum weighted degree from the graph
			map__node__deleted[c_highest_priority_node_id] = true;
			//
			// update priorities on priority queue
			// iterate over all not removed neighbours of
			// c_highest_priority_node.
			list_of_all_original_neighbours_of_c_highest_priority_node_id = this.graph_as_adjacency_list[c_highest_priority_node_id];
			list_of_weights_of_all_original_neighbours_of_c_highest_priority_node_id = this.graph_as_adjacency_list_weights[c_highest_priority_node_id];
			for (i = 0; i < list_of_all_original_neighbours_of_c_highest_priority_node_id.length; i++) {
				c_not_removed_neighbour_id = list_of_all_original_neighbours_of_c_highest_priority_node_id[i];
				if (map__node__deleted[c_not_removed_neighbour_id]) {
					continue;
				}
				//
				// Update priority of the not_removed_neighbour_of_c_min_weighted_degree_node
				PRIORITY_of_c_not_removed_neighbour_of_c_min_weighted_degree_node = pq
						.keyOf(c_not_removed_neighbour_id);
				weight_on_current_edge = list_of_weights_of_all_original_neighbours_of_c_highest_priority_node_id[i];
				new_PRIORITY_of_c_not_removed_neighbour_of_c_min_weighted_degree_node = PRIORITY_of_c_not_removed_neighbour_of_c_min_weighted_degree_node
						- weight_on_current_edge;
				//

				map__node_id__load[c_highest_priority_node_id] += weight_on_current_edge;

				//
				pq.updateKey(c_not_removed_neighbour_id,
						new_PRIORITY_of_c_not_removed_neighbour_of_c_min_weighted_degree_node);
				//
				// update c_E_w_value
				c_E_w_value -= weight_on_current_edge;
			}
			//
			// update c_V_value
			c_V_value--;
		}
		//
		//
		this.output_best_V_size = best_solution__size;
		this.output_best_solution_value = best_solution__value;
		//
		return;
	}
}
/*
 * 
 * Further, researchers may be interested in obtaining patterns containing pre-
 * specified nodes. We refer to this as the subset maximum density problem, and
 * this is described in Section 3.
 * 
 * In [13], Saha et al. give an O(n4 log n)-time algorithm for the densest sub-
 * graph with a specified subset problem based on the maximum flow (min-cut)
 * problem. Later, Wenbin et al. propose a LP-based polynomial time algorithm
 * for this problem with less time O(n3.5) [3]. It is an open problem whether
 * there exist any less time exact algorithms or approximation algorithms for
 * the densest subgraph with a specified subset problem.
 * 
 * 
 * We also show that the densest subgraph problem with a specified subset of
 * vertices that have to be included in the solution can be solved optimally in
 * polynomial time.
 */

/*
 * 3.1.6. Densest Subgraph with Specified Subset Motivated by gene annotation
 * graphs in computational biology, Saha et al. [13] introduced the problem of
 * finding a maximum edge density subgraph, under the constraint that it
 * contains a specified subset of nodes. While this constraint makes the problem
 * harder, it remains solvable in polynomial time. Saha et al. [13] gave an O(n
 * 4 log n)-time algorithm for the problem in graphs with n vertices, using
 * network flows. Later, Chen et al. [14] found an O(n 2 )-time greedy
 * approximation algorithm with approximation ratio 2(1 + k 3 ), for the case
 * where the subgraph is required to contain at least k vertices. The results
 * carry over to the edge weighted case, as well.
 * 
 * 
 */
