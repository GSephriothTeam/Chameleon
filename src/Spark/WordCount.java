package Spark;

import org.apache.spark.api.java.*;
import org.apache.spark.SparkConf;

/**
 *
 * Created by xuzhuchen on 11/19/17.
 */
/* SimpleApp.java */
import org.apache.spark.api.java.*;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.function.Function;

public class WordCount {
    public static void main(String[] args) {

        String logFile = "/Data/Data10_1000.csv"; // Should be some file on your system
        SparkConf conf = new SparkConf().setAppName("CountData").setMaster("local");
        JavaSparkContext sc = new JavaSparkContext(conf);
        JavaRDD<String> logData = sc.textFile(logFile).cache();

        long numRows = logData.filter(new Function<String, Boolean>() {
            public Boolean call(String s) { return s.contains(","); }
        }).count();

        System.out.println("There are "+numRows+" rows.");

        sc.stop();
    }
}