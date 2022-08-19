const
	ε = 1e-6,
	δ = Math.abs,
	Y = 1.790857183,
	K1 = 1.0148,
	K2 = 0.23185,
	K3 = -0.14499,
	K4 = 0.02406,
	C1 = K1,
	C2 = 5 * K2,
	C3 = 7 * K3,
	C4 = 9 * K4

export function patterson(λ, φ) {
	return [λ, φ * (K1 + φ ** 4 * (K2 + φ ** 2 * (K3 + K4 * φ ** 2)))];
}

patterson.invert = function (x, y) {
	if (y > Y) y = Y;
	else if (y < -Y) y = -Y;
	var yc = y, delta;

	do { // Newton-Raphson
		var y2 = yc * yc;
		yc -= delta = ((yc * (K1 + y2 * y2 * (K2 + y2 * (K3 + K4 * y2)))) - y) / (C1 + y2 * y2 * (C2 + y2 * (C3 + C4 * y2)));
	} while (δ(delta) > ε);

	return [x, yc];
}
