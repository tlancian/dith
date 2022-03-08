package tools;

import java.util.NoSuchElementException;

public class FastAndSimpleIndexMinPQ {

	protected int maxN; // maximum number of elements on PQ
	protected int n; // number of elements on PQ
	protected int[] pq; // binary heap using 1-based indexing
	protected int[] qp; // inverse of pq - qp[pq[i]] = pq[qp[i]] = i
	protected double[] keys; // keys[i] = priority of i
	//
	protected int swap;
	//

	public FastAndSimpleIndexMinPQ(int maxN) {
		this.maxN = maxN;
		n = 0;
		keys = new double[maxN + 1]; // make this of length maxN??
		pq = new int[maxN + 1];
		qp = new int[maxN + 1]; // make this of length maxN??
		for (int i = 0; i <= maxN; i++) {
			qp[i] = -1;
		}
		swap = -1;
		return;
	}

	public void insert(int i, double key) {
		if (qp[i] != -1) {
			return;
		}
		n++;
		qp[i] = n;
		pq[n] = i;
		keys[i] = key;
		swim(n);
	}

	public double minKey() {
		if (n == 0)
			throw new NoSuchElementException("Priority queue underflow");
		return keys[pq[1]];
	}

	public int delMin() {
		if (n == 0)
			throw new NoSuchElementException("Priority queue underflow");
		int min = pq[1];
		// exch(1, n--);
		//
		swap = pq[1];
		pq[1] = pq[n];
		pq[n] = swap;
		qp[pq[1]] = 1;
		qp[pq[n]] = n;
		//
		n--;
		//
		sink(1);
		assert min == pq[n + 1];
		qp[min] = -1; // delete
		keys[min] = Double.POSITIVE_INFINITY;
		pq[n + 1] = -1; // not needed
		return min;
	}

	public double keyOf(int i) {
		return keys[i];
	}

	public void updateKey(int i, double key) {
		//
		if (keys[i] > key) {
			//
			// decrease key
			keys[i] = key;
			swim(qp[i]);
		} else if (keys[i] < key) {
			//
			// increase key
			keys[i] = key;
			sink(qp[i]);
		}
		//
		return;
	}

	/***************************************************************************
	 * General helper functions.
	 ***************************************************************************/
//	private boolean greater(int i, int j) {
//		// return keys[pq[i]].compareTo(keys[pq[j]]) > 0;
//		return keys[pq[i]] > keys[pq[j]];
//	}

//	private void exch(int i, int j) {
//		swap = pq[i];
//		pq[i] = pq[j];
//		pq[j] = swap;
//		qp[pq[i]] = i;
//		qp[pq[j]] = j;
//		return;
//	}

	/***************************************************************************
	 * Heap helper functions.
	 ***************************************************************************/
	private void swim(int k) {
//		while (k > 1 && greater(k / 2, k)) {
		while (k > 1 && (keys[pq[k / 2]] > keys[pq[k]])) {
//			exch(k, k / 2);
			// exch
			swap = pq[k];
			pq[k] = pq[k / 2];
			pq[k / 2] = swap;
			qp[pq[k]] = k;
			qp[pq[k / 2]] = k / 2;
			//
			k = k / 2;
		}
		return;
	}

	private void sink(int k) {
		while (2 * k <= n) {
			int j = 2 * k;
//			if (j < n && greater(j, j + 1)) {
			if (j < n && (keys[pq[j]] > keys[pq[j + 1]])) {
				j++;
			}
			// if (!greater(k, j)) {
			if (keys[pq[k]] <= keys[pq[j]]) {
				break;
			}
			//
			// exch(k, j);
			//
			swap = pq[k];
			pq[k] = pq[j];
			pq[j] = swap;
			qp[pq[k]] = k;
			qp[pq[j]] = j;
			//
			k = j;
		}
		return;
	}
}
