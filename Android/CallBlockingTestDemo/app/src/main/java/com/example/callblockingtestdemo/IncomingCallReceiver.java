package com.example.callblockingtestdemo;

import android.app.Notification;
import android.app.NotificationChannel;
import android.app.NotificationManager;
import android.app.PendingIntent;
import android.content.BroadcastReceiver;
import android.content.Context;
import android.content.Intent;
import android.os.Build;
import android.telephony.TelephonyManager;
import android.view.Gravity;
import android.widget.Toast;

import androidx.core.app.NotificationCompat;

import java.lang.reflect.Method;
import com.android.internal.telephony.ITelephony;

public class IncomingCallReceiver extends BroadcastReceiver {


    @Override
    public void onReceive(Context context, Intent intent) {

        ITelephony telephonyService;
        try {
            String state = intent.getStringExtra(TelephonyManager.EXTRA_STATE);
            String number = intent.getExtras().getString(TelephonyManager.EXTRA_INCOMING_NUMBER);

            if(state.equalsIgnoreCase(TelephonyManager.EXTRA_STATE_RINGING)){
                TelephonyManager tm = (TelephonyManager) context.getSystemService(Context.TELEPHONY_SERVICE);
                try {
                    Method m = tm.getClass().getDeclaredMethod("getITelephony");

                    m.setAccessible(true);
                    telephonyService = (ITelephony) m.invoke(tm);

                } catch (Exception e) {
                    e.printStackTrace();
                }

                SendNotification( context, number );
                Toast.makeText(context, "Ring " + number, Toast.LENGTH_LONG).show();

            }
            //if(state.equalsIgnoreCase(TelephonyManager.EXTRA_STATE_OFFHOOK)){
            //    Toast.makeText(context, "Answered " + number, Toast.LENGTH_LONG).show();
            //}
            //if(state.equalsIgnoreCase(TelephonyManager.EXTRA_STATE_IDLE)){
            //    Toast.makeText(context, "Idle "+ number, Toast.LENGTH_LONG).show();
            //}
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    void SendNotification(Context context, String messageBody)
    {
        String CHANNEL_ID = "com.example.callblockingtestdemo.channel";

        // Create the NotificationChannel, but only on API 26+ because
        // the NotificationChannel class is new and not in the support library
        //if (Build.VERSION.SDK_INT >= Build.VERSION_CODES.O) {
            int importance = NotificationManager.IMPORTANCE_MAX;
            NotificationChannel channel = new NotificationChannel(CHANNEL_ID , "channel", importance);
            channel.setDescription("description");
            // Register the channel with the system; you can't change the importance
            // or other notification behaviors after this
            NotificationManager notificationManager = context.getSystemService(NotificationManager.class);
            notificationManager.createNotificationChannel(channel);
        //}


        Notification newMessageNotification = new Notification.Builder( context, CHANNEL_ID)
                .setContentTitle("Reverse caller lookup")
                .setContentText(messageBody)
                .setSmallIcon(android.R.drawable.ic_menu_info_details)
                .setAutoCancel(true)
                .build();

        notificationManager.notify(0, newMessageNotification);
    }
}