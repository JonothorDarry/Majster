#include <stdio.h>
#include <vector>
#include <math.h>
using namespace std;
using ll=long long;
const int D=51, DD=31;
const ll MAX=100000000000ll, C=10000000, G=MAX/C+1, CP=C+50, GP=G+50;

int lpf[CP], prime_divs[CP];
vector<int> all_primes;
int numbers_per_prime_amount[C][D];
ll dp[D][GP], pw2[D];

ll prime_counter(ll x){
	for (int k=G; k>=1; k--){
		ll real_value = x/k;
		dp[1][k] = real_value;
		int amount_analyzed = 0;

		for (int amount_primes = 2; amount_primes < D; amount_primes++){ //Add small primes
			dp[amount_primes][k] = 0;
			for (int ite=0; ite<all_primes.size(); ite++){ 
				ll j = all_primes[ite];
				if (j*j > real_value || pw2[amount_primes-1]*j > real_value){
					if (amount_primes == 2) amount_analyzed = ite;
				       	break;
				}

				ll current_value = j;
				int added_divisors = 1;
				while (current_value <= real_value && added_divisors <= amount_primes){
					int previous_amount_divisors = amount_primes - added_divisors;
					ll previous_max_value = real_value/current_value;

					if (real_value/current_value >= C) dp[amount_primes][k] += dp[previous_amount_divisors][k*current_value];
					else {
						int amount_lower_equal = numbers_per_prime_amount[previous_max_value][previous_amount_divisors];
						dp[amount_primes][k] += amount_lower_equal;
					}

					current_value *= j;
					added_divisors++;
				}
			}
		}

		for (ll j=2; j*j <= real_value; j++){ //Add huge primes
			if (k*j <= G) dp[prime_divs[j] + 1][k] += std::max(dp[1][k*j] - amount_analyzed, 0ll);
			else dp[prime_divs[j] + 1][k] += std::max(numbers_per_prime_amount[real_value/j][1] - amount_analyzed, 0); 
		}

		dp[0][k] = 1; //for one
		for (int amount_primes = 2; amount_primes < D; amount_primes++){
			dp[amount_primes][k] /= amount_primes;
			dp[1][k] -= dp[amount_primes][k];
		}
		dp[1][k] -= dp[0][k];
		//printf ("%d: pi(%lld) = %lld\n", k, real_value, dp[1][k]);
	}
	return dp[1][1];
}

void shared_part(ll maximal_x){
	pw2[0] = 1;
	for (int i=1; i<D; i++) pw2[i] = pw2[i-1]*2;

	for (int i=1; i<C; i++) lpf[i] = i;
	for (int i=2; i<C; i++){
		if (lpf[i] != i) continue;
		all_primes.push_back(i);

		for (int j=i; j<C; j+=i){
			if (lpf[j] == j) lpf[j] = i;
		}
	}

	numbers_per_prime_amount[1][0] = 1;
	for (int i=2; i<C; i++){
		for (int j=0; j<DD; j++) numbers_per_prime_amount[i][j] = numbers_per_prime_amount[i-1][j];
	       	prime_divs[i] = prime_divs[i/lpf[i]] + 1;
		numbers_per_prime_amount[i][prime_divs[i]]++;
	}
}


int main(){
	shared_part(MAX);
	int t; scanf ("%d", &t);
	while(t--){
		ll x; scanf ("%lld", &x);
		printf ("%lld\n", prime_counter(x));
	}
return 0;}
