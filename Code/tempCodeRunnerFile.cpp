
db yh2o(db temp, db press, db molality) // Checked ......
{
    db co2 = kc02(temp, press, molality);
    db h2o = kh2o(temp, press);
    return ((1 - (1 / co2)) / ((1 / h2o) - (1 / co2)));
}

db xco2(db temp, db press, db molality)
{
    db yh20 = yh2o(temp, press, molality);
    db kco2 = kc02(temp, press, molality);