#include <stdio.h>
#include <vector>
#include <cmath>
#include "common.hpp"
using ll=long long;
using namespace std;
const ll MAX=100000000001ll, SQRT=ll(ceil(sqrt(MAX)));
const int BCv=25, Cv=1<<BCv, D=10001;

struct segtree_vertex{
	int value;
	int left;
	int right;
	int index;
};

vector <vector <segtree_vertex> > val(Cv<<1);
int unexplained[Cv];
int pi_value[Cv];
void initialize_segment_tree(ll limit){
	for (int i=0; i<Cv; i++) val[Cv + i].reserve(2);

	val[Cv].push_back({0, -1, -1, 0});
	val[Cv + 1].push_back({1, -1, -1, 0});
	for (int i=2; i<Cv; i++) val[Cv + i].push_back({0, -1, -1, 0});

	for (int i=1; i<Cv; i++) unexplained[i] = i;
}

void update_segment_tree(int basis, ll limit){
	for (ll value=basis; value <= limit; value+=basis){
		while (unexplained[value] % basis == 0) unexplained[value] /= basis;
		if (unexplained[value] > 1) continue;
		if (val[Cv + value].back().value == 1) continue;
		val[value+Cv].push_back({1, -1, -1, basis});
	}
}

void after_all_updates(){
	for (ll x=Cv-1; x>0; x--){
		int nexte = (x<<1);
		int l_length = val[nexte].size(), r_length = val[nexte+1].size();

		val[x].reserve(max(l_length, r_length)+1);
		for (int l=0,r=0; l<l_length || r<r_length; ){
			int index=-1;

			if (l < l_length && r < r_length && val[nexte+1][r].index == val[nexte][l].index) index = val[nexte][l].index, l++, r++;
			else if (l == l_length || (r < r_length && val[nexte+1][r].index < val[nexte][l].index)) index = val[nexte+1][r].index, r++;
			else index = val[nexte][l].index, l++;

			val[x].push_back({val[nexte][l-1].value + val[nexte+1][r-1].value, l-1, r-1, index});
		}
	}
}

ll query_segment_tree(ll l, ll r, int prime_nr){
	pair<int,int> node_left={1, prime_nr}, node_right = {1, prime_nr}; //(x,y)
	int res=0, bit = Cv>>1, l_root, r_root;

	int l_x, r_x;
	while (true){
		l_x = node_left.first, r_x = node_right.first;
		if (l_x%2 == 0 && r_x-l_x > 1) res += val[l_x+1][val[l_x>>1][l_root].right].value;
		if (r_x%2 == 1 && r_x-l_x > 1) res += val[r_x-1][val[r_x>>1][r_root].left].value;

		if (bit <= 0) break;
		l_root = node_left.second;
		r_root = node_right.second;

		if ((l&bit) == 0) node_left = {node_left.first<<1, val[node_left.first][node_left.second].left};
		else node_left = {(node_left.first<<1)+1, val[node_left.first][node_left.second].right};

		if ((r&bit) == 0) node_right = {node_right.first<<1, val[node_right.first][node_right.second].left};
		else node_right = {(node_right.first<<1)+1, val[node_right.first][node_right.second].right};

		bit>>=1;
	}
	res += val[l_x][node_left.second].value;
	if (l_x != r_x) res += val[r_x][node_right.second].value;

	return res;
}

void make_segment_tree(ll proper_limit, ll limit, vector <ll> &prime){
	initialize_segment_tree(proper_limit);
	for (int i=1; i<=limit; i++) update_segment_tree(prime[i], proper_limit);
	after_all_updates();
}

ll bases[D], tmp_bases[D];
void basis_of_the_psi(ll n, int width_of_antitree, int layers, vector <ll> &prime){
	for (ll j=1; j<=width_of_antitree; j++) bases[j] = 1;

	for (int i=1; i<=layers; i++){
		for (ll j=1; j <= width_of_antitree; j++){
			ll post_value = n/j;
			
			for (ll denominator = 1; denominator <= post_value; denominator*=prime[i]){
				if (j*denominator <= width_of_antitree) tmp_bases[j] += bases[j*denominator];
				else tmp_bases[j] += query_segment_tree(1, n/(j*denominator), i-1);
			}
		}
		for (int j=1; j<=width_of_antitree; j++) bases[j] = tmp_bases[j];
		for (int j=1; j<=width_of_antitree; j++) tmp_bases[j] = 0;
	}
}

ll calculate_psi(ll n, ll div, int a, int base_index, vector <ll> &prime){
	ll proper_res = bases[div];
	ll proper_value = n/div;

	for (int i=base_index+1; i<=a; i++){
		for (ll denominator = prime[i]; denominator <= proper_value; denominator*=prime[i]){
			proper_res += query_segment_tree(1, n/(div*denominator), i-1);
		}
	}
	return proper_res;
}


vector <ll> prime;
ll prime_counters[SQRT];
ll prime_counter(ll x){
	if (x==1) return 0;
	ll point_sqrt_3 = ll(ceil(pow(x, 1.0f/3.0f)));
	ll base_index = get_prime_above(prime, ll(ceil(pow(x, 1/3.)))); //to binary searcher
	basis_of_the_psi(x, point_sqrt_3, base_index, prime);
	
	for (int i=point_sqrt_3; i>0; i--){
		ll point = x/i;
		ll point_sqrt_2 = ll(ceil(sqrt(point)));
		ll alpha_index = max(base_index, get_prime_above(prime, point_sqrt_2));
		prime_counters[i] = point - calculate_psi(x, i, alpha_index, base_index, prime) + alpha_index;
		for (ll j=2; j<=point_sqrt_2; j++){
			if (j*i <= point_sqrt_3) prime_counters[i] -= max(prime_counters[i*j]-alpha_index, 0ll);
			else prime_counters[i] -= max(pi_value[point/j]-alpha_index, 0ll);
		}
		//printf ("B %lld %lld\n", point, prime_counters[i]);
	}

	return prime_counters[1];
}

void shared_part(ll maximal_x){
	ll last_prime = ll(ceil(pow(maximal_x, 2/3.)))+1;
	get_primes(last_prime, prime);
	ll limiting_prime = get_prime_above(prime, ll(ceil(pow(maximal_x, 1/3.)))); //to binary searcher
	ll limiting_prime_sqrt = get_prime_above(prime, ll(ceil(pow(maximal_x, 1/2.)))); //to binary searcher

	ll proper_limit = maximal_x/prime[limiting_prime];
	make_segment_tree(proper_limit, limiting_prime_sqrt, prime);
	for (auto &x: prime) pi_value[x]++;

	pi_value[0] = 0;
	for (int i=1; i<Cv; i++) pi_value[i] += pi_value[i-1];
}

int main(){
	shared_part(MAXIM);

	int t; scanf ("%d", &t);
	while(t--){
		ll x; scanf ("%lld", &x);
		printf ("%lld\n", prime_counter(x));
	}
return 0;}
