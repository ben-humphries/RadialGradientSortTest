#include <SFML/Graphics.hpp>
#include <math.h>
#include <thread>
#include <chrono>
#include <algorithm>

#define PI 3.14159256

typedef struct {
	double r;       // a fraction between 0 and 1
	double g;       // a fraction between 0 and 1
	double b;       // a fraction between 0 and 1
} rgb;
typedef struct {
	double h;       // angle in degrees
	double s;       // a fraction between 0 and 1
	double v;       // a fraction between 0 and 1
} hsv;
static rgb   hsv2rgb(hsv in);

void pollEvents();
void updateCircle();
void updateColors();

void shuffle();
void wait(long milliseconds);
void swap(int * a, int * b);
void flip(int arr[], int i);
int findMax(int arr[], int n);
void merge(int l, int m, int r);

void selectionSort();
void bubbleSort();
void insertionSort();
void quickSort(int low, int high);
void shellSort();
void cycleSort();
void pancakeSort();
void cocktailSort();
void mergeSort(int low, int high);

clock_t t;
float elapsedTime = 0;

const int numTriangles = 3000;
int radius = 300;
const int delay = 0; //ms

float colorIncrement = 360.0 / (numTriangles);

float theta = 0;
float dtheta = (2 * PI) / numTriangles;

sf::Vector2f center = sf::Vector2f(300, 300);

sf::RenderWindow window(sf::VideoMode(600, 600), "Selection Sort");

sf::Vertex vertices[numTriangles * 2 + 1];
sf::Color colors[numTriangles];
int sortArray[numTriangles];

bool loop = true;

int main()
{

	srand(time(0));

	for (int i = 0; i < numTriangles; i++) {
		sortArray[i] = i;
	}

	shuffle();

	updateColors();

	vertices[0] = sf::Vertex(center, sf::Color::White);

	updateCircle();

	while (window.isOpen())
	{
		//Calculate dt
		float dt = ((float)(clock() - t)) / CLOCKS_PER_SEC;
		t = clock();
		elapsedTime += dt;

		if (elapsedTime > 2 && loop) {
			//selectionSort();
			//insertionSort();
			//bubbleSort();;
			//quickSort(0, numTriangles-1);
			shellSort();
			//cycleSort();
			//pancakeSort();
			//cocktailSort();
			//mergeSort(0, numTriangles-1);
			updateColors();
			updateCircle();

			loop = false;
		}

		pollEvents();

		window.clear(); 
		window.draw(vertices, numTriangles * 2 + 1, sf::TrianglesFan);
		window.display();
	}

	return 0;
}

void pollEvents() {
	sf::Event event;
	while (window.pollEvent(event))
	{
		if (event.type == sf::Event::Closed)
			window.close();
	}
}

void updateCircle() {
	theta = 0;

	for (int i = 1; i < numTriangles * 2 + 1; i += 2) {

		vertices[i] = sf::Vertex(center + sf::Vector2f(radius*cos(theta), radius*sin(theta)), colors[sortArray[i / 2]]);
		vertices[i + 1] = sf::Vertex(center + sf::Vector2f(radius*cos(theta + dtheta), radius*sin(theta + dtheta)), colors[sortArray[i / 2]]);

		theta += dtheta;
	}
}

void updateColors() {

	float hue = 0;
	for (int i = 0; i < numTriangles; i++) {

		hsv a;
		a.h = hue;
		a.s = 1;
		a.v = 1;

		rgb b = hsv2rgb(a);

		colors[i] = sf::Color(b.r * 255, b.g * 255, b.b * 255);

		hue += colorIncrement;
	}
}
void shuffle() {
	std::random_shuffle(std::begin(sortArray), std::end(sortArray));
}

void wait(long milliseconds) {
	std::this_thread::sleep_for(std::chrono::milliseconds(milliseconds));
}

void swap(int * a, int * b) {
	int t = *a;
	*a = *b;
	*b = t;

	colors[*a] = sf::Color::White;

	updateCircle();
	updateColors();
	window.clear();
	window.draw(vertices, numTriangles * 2 + 1, sf::TrianglesFan);
	window.display();
	wait(delay);
}


void selectionSort() {

	int i, j, min_idx;

	for (i = 0; i < numTriangles - 1; i++)
	{
		min_idx = i;
		for (j = i + 1; j < numTriangles; j++)
			if (sortArray[j] < sortArray[min_idx])
				min_idx = j;
		swap(&sortArray[min_idx], &sortArray[i]);

		pollEvents();
	}

}

void bubbleSort()
{
	int i, j;
	for (i = 0; i < numTriangles - 1; i++)

		for (j = 0; j < numTriangles - i - 1; j++)
			if (sortArray[j] > sortArray[j + 1]) {
				swap(&sortArray[j], &sortArray[j+1]);

				pollEvents();
			}
}

void insertionSort()
{
	int i, key, j;
	for (i = 1; i < numTriangles; i++)
	{
		key = sortArray[i];
		j = i - 1;

		while (j >= 0 && sortArray[j] > key)
		{
			sortArray[j + 1] = sortArray[j];
			j = j - 1;
		}
		sortArray[j + 1] = key;

		colors[sortArray[i]] = sf::Color::White;

		updateCircle();
		updateColors();
		window.clear();
		window.draw(vertices, numTriangles * 2 + 1, sf::TrianglesFan);
		window.display();
		wait(delay);

		pollEvents();
	}
}


void quickSort(int low, int high)
{
	pollEvents();

	if (low < high)
	{
		/* pi is partitioning index, arr[p] is now
		at right place */
		int pivot = sortArray[high];    // pivot
		int i = (low - 1);  // Index of smaller element

		for (int j = low; j <= high - 1; j++)
		{
			// If current element is smaller than or
			// equal to pivot
			if (sortArray[j] <= pivot)
			{
				i++;    // increment index of smaller element
				swap(&sortArray[i], &sortArray[j]);

			}
		}
		swap(&sortArray[i+1], &sortArray[high]);

		int pi = (i + 1);

		// Separately sort elements before
		// partition and after partition
		quickSort(low, pi - 1);
		quickSort(pi + 1, high);
	}
}

void shellSort()
{
	// Start with a big gap, then reduce the gap
	for (int gap = numTriangles / 2; gap > 0; gap /= 2)
	{
		// Do a gapped insertion sort for this gap size.
		// The first gap elements a[0..gap-1] are already in gapped order
		// keep adding one more element until the entire array is
		// gap sorted 
		for (int i = gap; i < numTriangles; i += 1)
		{
			// add a[i] to the elements that have been gap sorted
			// save a[i] in temp and make a hole at position i
			int temp = sortArray[i];

			// shift earlier gap-sorted elements up until the correct 
			// location for a[i] is found
			int j;
			for (j = i; j >= gap && sortArray[j - gap] > temp; j -= gap)
				sortArray[j] = sortArray[j - gap];

			//  put temp (the original a[i]) in its correct location
			sortArray[j] = temp;

			colors[sortArray[i]] = sf::Color::White;

			updateCircle();
			updateColors();
			window.clear();
			window.draw(vertices, numTriangles * 2 + 1, sf::TrianglesFan);
			window.display();
			wait(delay);

			pollEvents();
		}
	}
}

void cycleSort()
{

	for (int cycle_start = 0; cycle_start <= numTriangles - 2; cycle_start++)
	{
		int item = sortArray[cycle_start];

		int pos = cycle_start;
		for (int i = cycle_start + 1; i<numTriangles; i++)
			if (sortArray[i] < item)
				pos++;

		// If item is already in correct position
		if (pos == cycle_start)
			continue;

		// ignore all duplicate  elements
		while (item == sortArray[pos])
			pos += 1;

		// put the item to it's right position
		if (pos != cycle_start)
		{
			swap(&item, &sortArray[pos]);
		}

		// Rotate rest of the cycle
		while (pos != cycle_start)
		{
			pollEvents();

			pos = cycle_start;

			// Find position where we put the element
			for (int i = cycle_start + 1; i<numTriangles; i++)
				if (sortArray[i] < item)
					pos += 1;

			// ignore all duplicate  elements
			while (item == sortArray[pos])
				pos += 1;

			// put the item to it's right position
			if (item != sortArray[pos])
			{
				swap(&item, &sortArray[pos]);
			}
		}
	}

}

void pancakeSort()
{
	for (int curr_size = numTriangles; curr_size > 1; --curr_size)
	{

		int mi = findMax(sortArray, curr_size);

		if (mi != curr_size - 1)
		{

			flip(sortArray, mi);
			flip(sortArray, curr_size - 1);
		}
	}
}

void flip(int arr[], int i)
{
	int start = 0;
	while (start < i)
	{
		swap(&arr[start], &arr[i]);
		start++;
		i--;
	}
}

int findMax(int arr[], int n)
{
	int mi, i;
	for (mi = 0, i = 0; i < n; ++i)
		if (arr[i] > arr[mi])
			mi = i;
	return mi;
}

void cocktailSort()
{
	bool swapped = true;
	int start = 0;
	int end = numTriangles - 1;

	while (swapped)
	{
		swapped = false;

		for (int i = start; i < end; ++i)
		{
			if (sortArray[i] > sortArray[i + 1])
			{
				swap(&sortArray[i], &sortArray[i + 1]);
				swapped = true;
			}
		}
		if (!swapped)
			break;
		swapped = false;
		--end;

		for (int i = end - 1; i >= start; --i)
		{
			if (sortArray[i] > sortArray[i + 1])
			{
				swap(&sortArray[i], &sortArray[i + 1]);
				swapped = true;
			}
		}
		++start;
	}
}

void merge(int l, int m, int r)
{
	int i, j, k;
	int n1 = m - l + 1;
	int n2 = r - m;

	int * L = new int[n1];
	int * R = new int[n2];

	for (i = 0; i < n1; i++)
		L[i] = sortArray[l + i];
	for (j = 0; j < n2; j++)
		R[j] = sortArray[m + 1 + j];

	i = 0; 
	j = 0;
	k = l; 
	while (i < n1 && j < n2)
	{
		if (L[i] <= R[j])
		{
			sortArray[k] = L[i];
			i++;
		}
		else
		{
			sortArray[k] = R[j];
			j++;
		}

		colors[sortArray[k]] = sf::Color::White;

		updateCircle();
		updateColors();
		window.clear();
		window.draw(vertices, numTriangles * 2 + 1, sf::TrianglesFan);
		window.display();
		wait(delay);

		k++;
	}

	while (i < n1)
	{
		sortArray[k] = L[i];
		i++;
		k++;

		colors[sortArray[k-1]] = sf::Color::White;

		updateCircle();
		updateColors();
		window.clear();
		window.draw(vertices, numTriangles * 2 + 1, sf::TrianglesFan);
		window.display();
		wait(delay);
	}

	while (j < n2)
	{
		sortArray[k] = R[j];
		j++;
		k++;

		colors[sortArray[k-1]] = sf::Color::White;

		updateCircle();
		updateColors();
		window.clear();
		window.draw(vertices, numTriangles * 2 + 1, sf::TrianglesFan);
		window.display();
		wait(delay);
	}

	delete L;
	delete R;
}


void mergeSort(int low, int high)
{
	if (low < high)
	{

		int m = low + (high - low) / 2;

		colors[sortArray[m]] = sf::Color::White;

		updateCircle();
		updateColors();
		window.clear();
		window.draw(vertices, numTriangles * 2 + 1, sf::TrianglesFan);
		window.display();
		wait(delay);

		mergeSort(low, m);
		mergeSort(m + 1, high);

		merge(low, m, high);
	}
}


rgb hsv2rgb(hsv in)
{
	double      hh, p, q, t, ff;
	long        i;
	rgb         out;

	if (in.s <= 0.0) {       // < is bogus, just shuts up warnings
		out.r = in.v;
		out.g = in.v;
		out.b = in.v;
		return out;
	}
	hh = in.h;
	if (hh >= 360.0) hh = 0.0;
	hh /= 60.0;
	i = (long)hh;
	ff = hh - i;
	p = in.v * (1.0 - in.s);
	q = in.v * (1.0 - (in.s * ff));
	t = in.v * (1.0 - (in.s * (1.0 - ff)));

	switch (i) {
	case 0:
		out.r = in.v;
		out.g = t;
		out.b = p;
		break;
	case 1:
		out.r = q;
		out.g = in.v;
		out.b = p;
		break;
	case 2:
		out.r = p;
		out.g = in.v;
		out.b = t;
		break;

	case 3:
		out.r = p;
		out.g = q;
		out.b = in.v;
		break;
	case 4:
		out.r = t;
		out.g = p;
		out.b = in.v;
		break;
	case 5:
	default:
		out.r = in.v;
		out.g = p;
		out.b = q;
		break;
	}
	return out;
}